"""Offline weather-window -> process-response -> yield evidence, not mediation fractions."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import netCDF4
import numpy as np


WEATHER = ("temp", "prec", "sw", "lw", "wind")
MONTH_DAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def model_date(anchor_year, day):
    """One-based day on the model's 365-day calendar, not Gregorian timedelta."""
    if int(day) < 1:
        raise ValueError("model day must be positive")
    year, doy = anchor_year + (int(day) - 1) // 365, (int(day) - 1) % 365 + 1
    for month, length in enumerate(MONTH_DAYS, 1):
        if doy <= length:
            return f"{year:04d}-{month:02d}-{doy:02d}"
        doy -= length
    raise ValueError("invalid model date")


def file_sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze_counterfactual(path):
    with netCDF4.Dataset(path, "r") as dataset:
        if dataset.getncattr("schema_version") != 2:
            raise ValueError("process-chain analysis requires schema 2 attribution outputs")
        attrs = {key: dataset.getncattr(key) for key in dataset.ncattrs()}
        data = {name: np.ma.asarray(var[:], dtype=float).filled(np.nan)
                for name, var in dataset.variables.items()}
        diagnostics = {name.removeprefix("factual_"): (var.units, var.diagnostic_kind)
                       for name, var in dataset.variables.items()
                       if name.startswith("factual_") and "diagnostic_kind" in var.ncattrs()
                       and name != "factual_harvest_event"}

    day = data["day"].astype(int)
    anchor = int(attrs["anchor_year"])
    mask = data["weather_replacement_window"] == 1
    active = data["ad_active_window"] == 1
    start, end = int(attrs["replacement_first_day"]), int(attrs["replacement_last_day"])
    if not np.array_equal(day, np.arange(1, len(day) + 1)) or not 1 <= start <= end <= len(day):
        raise ValueError("invalid replacement/day coordinates")
    if not np.array_equal(data["calendar_year"], anchor + (day - 1) // 365) or not np.array_equal(
            data["day_of_year"], (day - 1) % 365 + 1):
        raise ValueError("calendar coordinates disagree with the 365-day contract")
    if not np.array_equal(mask, (day >= start) & (day <= end)):
        raise ValueError("replacement mask does not match its declared window")
    if attrs["counterfactual_kind"] == "window_restore" and np.any(mask & ~active):
        raise ValueError("window restoration must stay within factual growth days")
    for prefix, key in (("factual", "factual_harvest_day"), ("reference", "counterfactual_harvest_day")):
        expected = [int(attrs[key])] if int(attrs[key]) > 0 else []
        if list(day[data[f"{prefix}_harvest_event"] == 1]) != expected:
            raise ValueError("daily harvest events disagree with harvest metadata")
    common = (data["factual_harvest_event"] == 0) & (data["reference_harvest_event"] == 0)
    if not (mask & active).any():
        raise ValueError("no factual growing days in the replacement window")
    if 0 < int(attrs["counterfactual_harvest_day"]) < start:
        raise ValueError("counterfactual harvest precedes the weather intervention")

    weather, projected = {}, 0.0
    for name in WEATHER:
        factual, reference, gradient = (data[f"factual_{name}"], data[f"reference_{name}"],
                                        data[f"dyield_d{name}"])
        if not all(np.isfinite(x).all() for x in (factual, reference, gradient)):
            raise ValueError(f"nonfinite weather/gradient: {name}")
        if not np.array_equal(factual[~mask], reference[~mask]):
            raise ValueError(f"weather changed outside the declared window: {name}")
        if np.any(gradient[~active] != 0):
            raise ValueError("gradient is nonzero outside the fixed-event AD contract")
        delta = reference - factual
        contribution = float(np.dot(gradient, delta))
        weather[name] = {
            "reference_minus_factual_mean": float(delta[mask & active].mean()),
            "linearized_yield_change": contribution,
        }
        if name == "prec":
            weather[name]["reference_minus_factual_total_mm"] = float(delta[mask & active].sum())
        projected += contribution

    processes = {}
    for name, (units, kind) in diagnostics.items():
        factual, reference = data[f"factual_{name}"], data[f"reference_{name}"]
        if not common.any():
            # An immediate counterfactual harvest/failure can leave no matched
            # growth day. Preserve the event and yield; do not invent zero
            # process effects or discard every other window in this case.
            processes[name] = dict(units=units, kind=kind, first_resolved_difference_date=None,
                first_resolved_difference_lag_days=None, peak_signed_change=None, peak_date=None,
                mean_change_common_growth=None, end_common_growth_change=None)
            if kind == "daily_flux":
                processes[name]["integrated_change_common_growth"] = None
            continue
        if not np.isfinite(factual[common]).all() or not np.isfinite(reference[common]).all():
            raise ValueError(f"nonfinite process value on common growing days: {name}")
        # Prefix differences indicate a broken paired experiment, not a lagged
        # physical response. Ignore NaN padding and harvest-reset values.
        prefix = common & (day < start)
        if not np.array_equal(factual[prefix], reference[prefix]):
            raise ValueError(f"process changed before the weather intervention: {name}")
        delta = reference - factual
        changed = common & ~np.isclose(reference, factual, rtol=1e-5, atol=1e-7)
        first = int(day[np.flatnonzero(changed)[0]]) if changed.any() else None
        indices = np.flatnonzero(common)
        peak = int(indices[np.argmax(np.abs(delta[indices]))])
        processes[name] = {
            "units": units, "kind": kind,
            "first_resolved_difference_date": model_date(anchor, first) if first is not None else None,
            "first_resolved_difference_lag_days": first - start if first is not None else None,
            "peak_signed_change": float(delta[peak]), "peak_date": model_date(anchor, day[peak]),
            "mean_change_common_growth": float(delta[common].mean()),
            "end_common_growth_change": float(delta[indices[-1]]),
        }
        if kind == "daily_flux":
            processes[name]["integrated_change_common_growth"] = float(delta[common].sum())

    factual_yield = float(attrs["factual_yield"])
    counterfactual_yield = float(attrs["counterfactual_yield"])
    harvest_day = int(attrs["counterfactual_harvest_day"])
    observed = harvest_day > 0
    if not np.isfinite(factual_yield) or (observed and not np.isfinite(counterfactual_yield)):
        raise ValueError("observed harvest has nonfinite yield")
    if not observed and not np.isnan(counterfactual_yield):
        raise ValueError("unharvested counterfactual must have missing yield, not zero")
    yield_change = counterfactual_yield - factual_yield if observed else None
    management_inputs = {}
    for name in ("fertilizer_input", "manure_input"):
        if name not in diagnostics:
            continue
        management_inputs[name] = {
            prefix: [{"date": model_date(anchor, day[index]), "gN_m2": float(value)}
                     for index, value in enumerate(data[f"{prefix}_{name}"])
                     if np.isfinite(value) and value > 0]
            for prefix in ("factual", "reference")
        }
    # Day-level ordering here is diagnostic. A first difference is a numerical
    # resolution threshold, not an inferred causal onset or a significance test.
    summary = {
        "file": Path(path).name, "cell_id": int(attrs["cell_id"]),
        "reference_year": int(attrs["reference_year"]), "kind": attrs["counterfactual_kind"],
        "replacement_start": model_date(anchor, start), "replacement_end": model_date(anchor, end),
        "weather_statistics_days": int((mask & active).sum()),
        "weather_statistics_scope": "replacement intersect factual growth days",
        "factual_fphu_range_in_window": [float(np.nanmin(data["factual_fphu"][mask & active])),
                                          float(np.nanmax(data["factual_fphu"][mask & active]))],
        "factual_yield": factual_yield, "counterfactual_yield": counterfactual_yield if observed else None,
        "yield_change": yield_change, "yield_units": attrs["yield_units"],
        "factual_harvest_date": model_date(anchor, attrs["factual_harvest_day"]),
        "counterfactual_harvest_date": model_date(anchor, harvest_day) if observed else None,
        "harvest_shift_days": harvest_day - int(attrs["factual_harvest_day"]) if observed else None,
        "crop_failed": bool(attrs["counterfactual_failed"]),
        "schedule_matches": bool(attrs["schedule_matches"]),
        "linearized_fixed_event_yield_change": projected,
        "linearization_residual": yield_change - projected if observed else None,
        "common_growth_days": int(common.sum()),
        "last_common_growth_date": model_date(anchor, day[common][-1]) if common.any() else None,
        "weather": weather, "processes": processes,
        "post_sowing_management_inputs": management_inputs,
    }
    return summary, data


def write_report(path, summaries, validation_status):
    def fmt(value):
        return "missing" if value is None else f"{value:+.4g}"

    lines = ["# Weather–window–process–yield evidence", "",
        "Positive changes mean reference-weather restoration minus factual weather.",
        "All five weather channels change jointly; per-variable AD terms are local approximations.",
        "Process comparisons use common growing days only, excluding harvest resets and NaN padding.",
        "PHU fraction is model progress, not an observed flowering stage.",
        "First resolved differences use rtol=1e-5 / atol=1e-7, not statistical significance.",
        "Single-window finite yield changes are NOT additive process contributions.",
        "Full-reference runs may also change weather after the factual harvest; window runs do not.",
        f"Case directional finite-difference status: `{validation_status}`.", "",
        "| Reference | Restored window | Kind | Yield change (t DM/ha) | Harvest shift (days) | Failed | ΔGPP common days (gC/m²) |",
        "| --- | --- | --- | ---: | ---: | --- | ---: |"]
    for item in summaries:
        lines.append(f"| {item['reference_year']} | {item['replacement_start']} → {item['replacement_end']} | "
            f"{item['kind']} | {fmt(item['yield_change'])} | {fmt(item['harvest_shift_days'])} | "
            f"{item['crop_failed']} | {fmt(item['processes']['gpp']['integrated_change_common_growth'])} |")
    for item in summaries:
        lines += ["", f"## {item['file']}", "",
            f"Weather restoration: ΔT = {fmt(item['weather']['temp']['reference_minus_factual_mean'])} °C; "
            f"Δprecipitation = {fmt(item['weather']['prec']['reference_minus_factual_total_mm'])} mm "
            "during the factual-growth part of the window.",
            f"Fixed-event first-order yield change = {fmt(item['linearized_fixed_event_yield_change'])}; "
            f"ordinary finite change minus this approximation = {fmt(item['linearization_residual'])} t DM/ha.",
            "A shifted harvest/failure is not covered by the fixed-event derivative.", "",
            "| Process diagnostic | First resolved difference | Lag (days) | Peak signed change | Units |",
            "| --- | --- | ---: | ---: | --- |"]
        for name in ("temperature_stress", "rootzone_available_water", "water_sufficiency",
                     "nitrogen_uptake", "nitrogen_limitation", "gpp", "npp", "fphu", "storage_carbon"):
            process = item["processes"][name]
            lines.append(f"| {name} | {process['first_resolved_difference_date'] or 'none resolved'} | "
                f"{fmt(process['first_resolved_difference_lag_days'])} | {fmt(process['peak_signed_change'])} | {process['units']} |")
        for name, inputs in item["post_sowing_management_inputs"].items():
            lines += ["", f"Post-sowing {name} (date, gN/m²): factual {inputs['factual']}; reference {inputs['reference']}."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_trace(path, summary, data):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    panels = (("temp", "Temperature (°C)"), ("prec", "Precipitation (mm/day)"),
              ("rootzone_available_water", "Root-weighted water\n(top 3 layers, mm)"),
              ("nitrogen_limitation", "N-supported Vcmax\nfraction"),
              ("gpp", "GPP (gC/m²/day)"), ("storage_carbon", "Storage C (gC/m²)"))
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 12), sharex=True, constrained_layout=True)
    for axis, (name, label) in zip(axes, panels):
        for prefix, color in (("factual", "#ad3939"), ("reference", "#187b74")):
            visible = (data[f"{prefix}_harvest_event"] == 0)
            # Weather remains visible on each observed harvest day; crop reset
            # states do not. Do not draw lines through missing/no-harvest data.
            if name in WEATHER:
                visible |= data[f"{prefix}_harvest_event"] == 1
            axis.plot(data["day"][visible], data[f"{prefix}_{name}"][visible], label=prefix, color=color, lw=1.4)
        window = data["day"][data["weather_replacement_window"] == 1]
        growing = data["day"][(data["factual_harvest_event"] == 0) | (data["reference_harvest_event"] == 0)]
        axis.axvspan(window[0] - 0.5, min(window[-1], growing[-1]) + 0.5, color="#cccccc", alpha=0.25)
        axis.set_xlim(growing[0] - 1, growing[-1] + 1)
        axis.set_ylabel(label)
        axis.grid(alpha=0.15)
    axes[0].legend(loc="best")
    axes[-1].set_xlabel("Joined model day (365-day calendar; dates in report)")
    yield_label = "missing" if summary["yield_change"] is None else f"{summary['yield_change']:+.4f}"
    fig.suptitle(f"{summary['kind']} · reference {summary['reference_year']} · "
                 f"Δyield={yield_label} t DM/ha · harvest shift={summary['harvest_shift_days']} days · failed={summary['crop_failed']}\n"
                 "Paired process responses, not additive causal process fractions", fontsize=11)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    parser.add_argument("--output", type=Path, help="new directory inside the isolated Paper 1 results root")
    parser.add_argument("--plots", action="store_true", help="also render one static diagnostic figure per counterfactual")
    args = parser.parse_args()
    private_environment = Path(__file__).resolve().parents[3] / "environment"
    if not Path(sys.prefix).resolve().is_relative_to(private_environment.resolve()):
        parser.error("use the dedicated attribution/environment Python environment")
    directory = args.case_directory.resolve()
    roots = [p for p in directory.parents if p.name == "default_parameter_attribution_v1"
             and p.parent.name == "output" and p.parent.parent.name == "paper1_global_cell"]
    if len(roots) != 1:
        parser.error("input must be inside the isolated Paper 1 results root")
    root = roots[0]
    output = (args.output or directory / "process_chain_analysis").resolve()
    if output == root or not output.is_relative_to(root) or output.exists():
        parser.error("output must be a new directory inside the isolated results root")
    manifest_path = directory / "study_complete.toml"
    with manifest_path.open("rb") as stream:
        manifest = tomllib.load(stream)
    if manifest.get("status") != "complete" or manifest.get("stage") != "attribution":
        parser.error("only completed attribution cases can be analyzed")
    records = []
    for reference in manifest["references"]:
        records.append(reference)
        records.extend(reference.get("windows", []))
    if not records or not any(r.get("counterfactual_kind") == "window_restore" for r in records):
        parser.error("no paired weather-window counterfactuals in this case")
    summaries, paths = [], []
    for record in records:
        # Server absolute paths are provenance. Downloaded case files retain
        # their basenames, and their hashes must still match the manifest.
        path = (directory / Path(record["output_path"]).name).resolve()
        if not path.is_relative_to(directory) or file_sha(path) != record["output_sha256"]:
            parser.error(f"counterfactual output is missing, escaped or changed: {path}")
        summary, _ = analyze_counterfactual(path)
        if path in paths or summary["reference_year"] != record["reference_year"] or summary["kind"] != record["counterfactual_kind"]:
            parser.error("duplicate or mismatched counterfactual record")
        if summary["cell_id"] != manifest["case"]["cell_id"] or not np.isclose(
                summary["factual_yield"], manifest["production_yield"], rtol=1e-6, atol=1e-7):
            parser.error("counterfactual does not belong to this factual case")
        summary["input_sha256"] = record["output_sha256"]
        summaries.append(summary)
        paths.append(path)
    output.mkdir(parents=True, exist_ok=False)
    report = {"scope": "paired weather-window process evidence; not additive mediation",
        "source_commit": manifest["repository_commit"], "manifest_sha256": file_sha(manifest_path),
        "analysis_script_sha256": file_sha(__file__), "python_version": sys.version,
        "numpy_version": np.__version__, "netcdf4_version": netCDF4.__version__,
        "directional_validation": manifest.get("directional_validation", {}), "counterfactuals": summaries}
    (output / "process_chain.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_report(output / "process_chain.md", summaries,
                 report["directional_validation"].get("status", "not recorded"))
    if args.plots:
        os.environ["MPLCONFIGDIR"] = str(private_environment / "matplotlib_cache")
        for path, summary in zip(paths, summaries):
            _, data = analyze_counterfactual(path)
            plot_trace(output / f"{path.stem}.png", summary, data)
    print(output)


if __name__ == "__main__":
    main()
