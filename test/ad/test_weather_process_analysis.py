"""Small synthetic NetCDF tests; no global input or parameter optimization."""

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import netCDF4
import numpy as np


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/attribution/analyze_weather_case.py"
SPEC = importlib.util.spec_from_file_location("weather_process_analysis", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def write_fixture(path):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("day", 730)
        ds.setncatts(dict(schema_version=2, anchor_year=2000, cell_id=30,
            reference_year=2001, counterfactual_kind="window_restore",
            replacement_first_day=361, replacement_last_day=368,
            factual_yield=2.0, counterfactual_yield=2.5, yield_units="t dry matter ha-1",
            factual_harvest_day=380, counterfactual_harvest_day=382,
            counterfactual_failed=0, schedule_matches=0))
        day = np.arange(1, 731)
        ds.createVariable("day", "i4", ("day",))[:] = day
        ds.createVariable("calendar_year", "i4", ("day",))[:] = 2000 + (day - 1) // 365
        ds.createVariable("day_of_year", "i4", ("day",))[:] = (day - 1) % 365 + 1
        mask = (day >= 361) & (day <= 368)
        active = (day >= 350) & (day < 380)
        ds.createVariable("weather_replacement_window", "i1", ("day",))[:] = mask
        ds.createVariable("ad_active_window", "i1", ("day",))[:] = active
        for name in analysis.WEATHER:
            ds.createVariable(f"factual_{name}", "f8", ("day",))[:] = 5.0
            ds.createVariable(f"reference_{name}", "f8", ("day",))[:] = 5.0 + mask * 0.5
            ds.createVariable(f"dyield_d{name}", "f8", ("day",))[:] = active * 0.2
        diagnostics = ("gpp", "npp", "rootzone_available_water", "nitrogen_limitation",
                       "temperature_stress", "water_sufficiency", "nitrogen_uptake", "fphu", "storage_carbon")
        for prefix, harvest in (("factual", 380), ("reference", 382)):
            growth = (day >= 350) & (day < harvest)
            event = np.full(730, np.nan)
            event[growth], event[harvest - 1] = 0, 1
            ds.createVariable(f"{prefix}_harvest_event", "f8", ("day",))[:] = event
            for name in diagnostics:
                values = np.full(730, np.nan)
                values[growth] = 5.0
                if prefix == "reference":
                    values[growth & (day >= 364)] += 0.2
                values[harvest - 1] = 0.0  # Must never become a false stress signal.
                var = ds.createVariable(f"{prefix}_{name}", "f8", ("day",))
                flux = name in ("gpp", "npp", "nitrogen_uptake")
                var.units = "gC m-2 day-1" if flux else "1"
                var.diagnostic_kind = "daily_flux" if flux else "state_or_auxiliary"
                var[:] = values


class ProcessChainTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.path = Path(self.scratch.name) / "fixture.nc"
        write_fixture(self.path)

    def test_noleap_dates(self):
        self.assertEqual(analysis.model_date(2000, 60), "2000-03-01")
        self.assertEqual(analysis.model_date(2000, 365), "2000-12-31")
        self.assertEqual(analysis.model_date(2000, 366), "2001-01-01")

    def test_window_lag_flux_and_harvest_reset(self):
        summary, _ = analysis.analyze_counterfactual(self.path)
        self.assertEqual(summary["replacement_start"], "2000-12-27")
        self.assertEqual(summary["replacement_end"], "2001-01-03")
        self.assertEqual(summary["yield_change"], 0.5)
        self.assertEqual(summary["harvest_shift_days"], 2)
        self.assertEqual(summary["common_growth_days"], 30)
        self.assertAlmostEqual(summary["linearized_fixed_event_yield_change"], 4.0)
        self.assertAlmostEqual(summary["linearization_residual"], -3.5)
        gpp = summary["processes"]["gpp"]
        self.assertEqual(gpp["first_resolved_difference_date"], "2000-12-30")
        self.assertEqual(gpp["first_resolved_difference_lag_days"], 3)
        self.assertAlmostEqual(gpp["integrated_change_common_growth"], 0.2 * 16)
        storage = summary["processes"]["storage_carbon"]
        self.assertAlmostEqual(storage["peak_signed_change"], 0.2)
        self.assertNotIn("integrated_change_common_growth", storage)
        json.dumps(summary, allow_nan=False)
        report = Path(self.scratch.name) / "report.md"
        analysis.write_report(report, [summary], "synthetic_not_validated")
        self.assertIn("NOT additive", report.read_text())

    def test_noop_has_no_process_difference(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            for name in list(ds.variables):
                if name.startswith("reference_"):
                    ds[name][:] = ds[name.replace("reference_", "factual_", 1)][:]
            ds.counterfactual_yield, ds.counterfactual_harvest_day, ds.schedule_matches = 2.0, 380, 1
        summary, _ = analysis.analyze_counterfactual(self.path)
        self.assertEqual(summary["yield_change"], 0)
        self.assertEqual(summary["linearized_fixed_event_yield_change"], 0)
        self.assertTrue(all(p["first_resolved_difference_date"] is None for p in summary["processes"].values()))

    def test_no_harvest_stays_missing(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds.counterfactual_yield, ds.counterfactual_harvest_day = np.nan, -1
            ds["reference_harvest_event"][381] = 0
        summary, _ = analysis.analyze_counterfactual(self.path)
        self.assertIsNone(summary["yield_change"])
        self.assertIsNone(summary["linearization_residual"])
        json.dumps(summary, allow_nan=False)

    def test_immediate_failure_keeps_yield_and_missing_process_comparison(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds.replacement_first_day = 350
            ds["weather_replacement_window"][:] = (np.arange(1, 731) >= 350) & (np.arange(1, 731) <= 368)
            ds.counterfactual_yield, ds.counterfactual_harvest_day, ds.counterfactual_failed = 0.0, 350, 1
            ds["reference_harvest_event"][:] = np.nan
            ds["reference_harvest_event"][349] = 1
        summary, _ = analysis.analyze_counterfactual(self.path)
        self.assertTrue(summary["crop_failed"])
        self.assertEqual(summary["yield_change"], -2.0)
        self.assertEqual(summary["common_growth_days"], 0)
        self.assertIsNone(summary["processes"]["gpp"]["integrated_change_common_growth"])
        self.assertIsNone(summary["last_common_growth_date"])
        json.dumps(summary, allow_nan=False)
        analysis.write_report(Path(self.scratch.name) / "failed.md", [summary], "not_validated")

    def test_outside_window_weather_change_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds["reference_temp"][350] += 0.1
        with self.assertRaisesRegex(ValueError, "outside the declared window"):
            analysis.analyze_counterfactual(self.path)

    def test_prefix_process_change_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds["reference_gpp"][350] += 0.1
        with self.assertRaisesRegex(ValueError, "before the weather intervention"):
            analysis.analyze_counterfactual(self.path)

    def test_nonfinite_active_process_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds["reference_gpp"][365] = np.nan
        with self.assertRaisesRegex(ValueError, "nonfinite process"):
            analysis.analyze_counterfactual(self.path)

    def test_mask_mismatch_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds["weather_replacement_window"][365] = 0
        with self.assertRaisesRegex(ValueError, "replacement mask"):
            analysis.analyze_counterfactual(self.path)

    def test_capacity_not_integrated_as_flux(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            for prefix in ("factual", "reference"):
                var = ds.createVariable(f"{prefix}_vcmax", "f8", ("day",))
                var.units, var.diagnostic_kind = "gC m-2 day-1", "state_or_auxiliary"
                var[:] = ds[f"{prefix}_gpp"][:]
        summary, _ = analysis.analyze_counterfactual(self.path)
        self.assertNotIn("integrated_change_common_growth", summary["processes"]["vcmax"])

    def test_management_pulse_shift_is_visible(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            for prefix, pulse in (("factual", 364), ("reference", 365)):
                var = ds.createVariable(f"{prefix}_fertilizer_input", "f8", ("day",))
                var.units, var.diagnostic_kind = "gN m-2 day-1", "daily_flux"
                var[:] = 0.0
                var[pulse - 1] = 4.0
        summary, _ = analysis.analyze_counterfactual(self.path)
        pulses = summary["post_sowing_management_inputs"]["fertilizer_input"]
        self.assertEqual(pulses["factual"], [{"date": "2000-12-30", "gN_m2": 4.0}])
        self.assertEqual(pulses["reference"], [{"date": "2000-12-31", "gN_m2": 4.0}])
        self.assertEqual(summary["processes"]["fertilizer_input"]["integrated_change_common_growth"], 0.0)

    def test_harvest_metadata_mismatch_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds.counterfactual_harvest_day = 381
        with self.assertRaisesRegex(ValueError, "harvest metadata"):
            analysis.analyze_counterfactual(self.path)

    def test_calendar_mismatch_rejected(self):
        with netCDF4.Dataset(self.path, "a") as ds:
            ds["day_of_year"][365] = 366
        with self.assertRaisesRegex(ValueError, "365-day contract"):
            analysis.analyze_counterfactual(self.path)

    def test_cli_provenance_isolation_and_no_overwrite(self):
        root = Path(self.scratch.name) / "paper1_global_cell/output/default_parameter_attribution_v1"
        case = root / "attribution/batches/synthetic_case"
        case.mkdir(parents=True)
        window = case / "window.nc"
        full = case / "reference_2001.nc"
        shutil.copyfile(self.path, window)
        shutil.copyfile(self.path, full)
        with netCDF4.Dataset(full, "a") as ds:
            ds.counterfactual_kind = "full_reference"
            ds.replacement_first_day, ds.replacement_last_day = 350, 730
            ds["weather_replacement_window"][:] = np.arange(1, 731) >= 350
        (case / "study_complete.toml").write_text(f'''status = "complete"
stage = "attribution"
repository_commit = "synthetic_fixture_not_scientific"
production_yield = 2.0
[case]
cell_id = 30
[directional_validation]
status = "synthetic_not_validated"
[[references]]
reference_year = 2001
counterfactual_kind = "full_reference"
output_path = "/server/reference_2001.nc"
output_sha256 = "{analysis.file_sha(full)}"
[[references.windows]]
reference_year = 2001
counterfactual_kind = "window_restore"
output_path = "/server/window.nc"
output_sha256 = "{analysis.file_sha(window)}"
''')
        command = [sys.executable, str(SCRIPT), str(case)]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        output = case / "process_chain_analysis/process_chain.json"
        report = json.loads(output.read_text())
        self.assertEqual(len(report["counterfactuals"]), 2)
        self.assertEqual(report["counterfactuals"][1]["input_sha256"], analysis.file_sha(window))
        before = analysis.file_sha(output)
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        self.assertEqual(analysis.file_sha(output), before)
        outside = Path(self.scratch.name) / "not_in_study"
        self.assertNotEqual(subprocess.run(command + ["--output", str(outside)], capture_output=True).returncode, 0)
        self.assertFalse(outside.exists())
        with netCDF4.Dataset(window, "a") as ds:
            ds.counterfactual_yield = 99.0
        fresh = case / "new_analysis"
        result = subprocess.run(command + ["--output", str(fresh)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("changed", result.stderr)
        self.assertFalse(fresh.exists())

    @unittest.skipUnless(os.getenv("WEATHER_CHAIN_JULIA_FIXTURE"), "optional Julia-to-Python integration fixture")
    def test_julia_netcdf_roundtrip(self):
        summary, data = analysis.analyze_counterfactual(Path(os.environ["WEATHER_CHAIN_JULIA_FIXTURE"]))
        self.assertEqual(summary["kind"], "window_restore")
        self.assertTrue(np.isfinite(summary["yield_change"]))
        self.assertIn("mineralization", summary["processes"])
        self.assertNotIn("integrated_change_common_growth", summary["processes"]["vcmax"])
        if os.getenv("WEATHER_CHAIN_PLOT_PATH"):
            analysis.plot_trace(Path(os.environ["WEATHER_CHAIN_PLOT_PATH"]), summary, data)


if __name__ == "__main__":
    unittest.main()
