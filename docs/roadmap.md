# Agrocosm.jl implementation roadmap

This document records implementation-level acceptance criteria. The shorter
public roadmap is in `docs/src/development/roadmap.md`.

## 1. Completed scientific and numerical baseline

- LPJmL-informed C3/C4 photosynthesis, water-limited `lambda`, respiration,
  allocation, phenology, sowing, harvest, failed-crop termination, and residue
  routing.
- Five-layer soil water and temperature, snow, phase-change energy, surface
  litter, coupled soil C/N decomposition, mineral-N transformations, gaseous
  losses, and leaching.
- Audited daily process order and C/N/water/energy balance diagnostics.
- CPU/GPU process kernels and `Float32`/`Float64` support.
- Canonical lifecycle state (`prognostic`, `fluxes`, `auxiliary`, `inputs`,
  `events`, `workspace`, and `output`) separated from `ProcessModules`.
- Backend-independent checkpoints, high-level simulation API, one-day
  `transition_day!`, streamed selected output, memory estimation, and runtime
  benchmark.
- Finite agricultural warm-up that leaves production time, output, and balance
  ledgers untouched while retaining warmed state; eager and streamed forcing
  paths are numerically equivalent.

This foundation remains the scientific regression baseline. Alternative
processes must demonstrate their differences against it rather than silently
replacing it.

## 2. AgrocosmData status

Milestones 1–5 are substantially complete at the code and fixture-test level:

- dataset catalog and versioned backend-neutral contracts;
- canonical grid selection and compact/global round trips;
- 12-PFT registry, explicit 64/32/24/16-band mappings, and crop masks;
- soil-code properties, pH, sowing date, PHU, fertilizer, manure, residue, and
  land-use readers;
- HWSD SOC/total-N aggregation, vertical remapping, uncertainty/fallback
  provenance, field-capacity water, and native initial state;
- daily temperature, precipitation, net longwave, and downward shortwave
  streaming; annual CO₂ alignment; 365-day normalization; block prefetch;
- full ten-cell equivalence through `model_initial_data` and
  `climate_forcings`.
- configuration-driven extraction of the 2015 rainfed-wheat management fields
  and 2015–2016 daily climate for a bounded global test dataset.

Remaining data-layer work is production hardening rather than new loader
architecture:

- execute and quality-control the full canonical-grid HWSD product;
- preserve full source/provenance manifests for server runs;
- benchmark real server NetCDF access and add a canonical cache only if direct
  compact reads are too slow.

Warm-up, backend transfer, state evolution, and global
execution remain responsibilities of Agrocosm.jl, not AgrocosmData.jl.

## 3. Immediate production sequence

### 3.1 Fixed 2015 wheat domain

- Select compact cells where 2015 rainfed-wheat `landfrac > 0`.
- Use `landfrac` only for selection and provenance; do not multiply any model
  process or reported crop quantity by fractional area.
- Reuse the 2015 sowing date, PHU, fertilizer, manure, and residue settings in
  every simulated production year, matching the chosen ISIMIP experiment.

Acceptance: every selected cell runs one rainfed-wheat stand in stable compact
ordering, and changing a positive land fraction without changing its sign does
not alter a cell-level model trajectory.

### 3.2 Streamed agricultural warm-up

- Accept restartable complete-year climate-block readers without materializing
  a global year.
- Cycle one or more historical years for the configured warm-up duration.
- Preserve production `simulated_days == 0`, empty production output, and
  untouched production balance ledgers.
- Record annual litter/fast/slow/total C/N, mineral N, and water summaries and
  save the final native checkpoint.

Acceptance: streamed and eager warm-up are numerically equal, including the
final prognostic state and annual report.

Status: the streamed/eager implementation and equivalence regression are
complete. The production runner writes and exactly restores the warm-up state,
then repeats checkpoint/restart at the end of 2015.

### 3.3 Global rainfed-wheat smoke test

- The local subset contains fixed 2015 management and 2015–2016 forcing. Its
  first ten `landfrac > 0` cells pass a 730-day CPU smoke test with HWSD state.
- Generate the complete canonical-grid HWSD product and review its
  coverage and conservation summaries.
- Start with one rainfed wheat PFT over all cells selected by land use.
- Run CPU and a single GPU using identical compact cell ordering.
- Stream climate and monthly/annual output; avoid full daily global ledgers.
- Check NaN/Inf, invalid negative pools, crop lifecycle failures, memory peak,
  throughput, restart continuity, and sampled or online balance closure.
- Reconstruct outputs to `720 × 280` by `cellid` and verify mask alignment.
- Continue through the second forcing year to test cross-year state and
  checkpoint/restart continuity after the first-year smoke test passes.

Acceptance: the complete year finishes within estimated memory, CPU/GPU
differences meet declared tolerances, and the second-year restart/reassembly is
deterministic.

### 3.4 HWSD pool-allocation decision

The current native initialization conserves HWSD layer SOC and total N using a
documented 40:60 fast/slow split and zero litter. Do not add an elaborate
equilibrium allocator without evidence. First inspect the ten-year warm-up for
initial respiration pulses, litter/fast-pool stabilization, mineral-N drift,
and total C/N trajectories. If needed, implement a constrained allocation that
preserves every layer total and records uncertainty.

Current decision: retain 40:60 as the first global production baseline, not as
an equilibrium claim. In the real ten-cell ten-year warm-up, total C fell
13.3%, total N fell 11.1%, and the fast-C fraction moved from 0.400 to 0.314.
Year 10 still lost 1.27% C and 0.97% N, so the report correctly returns
`review_pool_allocation`. Because the total pools are still drifting, changing
only the initial fast/slow ratio is not a defensible fix; reconsider it together
with a longer or target-constrained spin-up after the global baseline run.

## 4. Differentiable daily transition

- Keep warm-up, I/O, diagnostics, and reporting outside the active path.
- Select a small continuous parameter/state set for the first Enzyme CPU
  smoke test.
- Compare gradients with finite differences on smooth, event-free windows.
- Classify bisection, min/max clamps, sowing, harvest, fertilization, and crop
  failure as smooth, piecewise, inactive, or requiring an alternative mode.
- Add GPU AD only after CPU primal and gradient regressions are stable.

## 5. Alternative canopy exchange

The current LPJmL-informed canopy remains the default global daily pathway.
Phase 3 may add a separately configured Farquhar–Medlyn–Penman–Monteith
alternative. It requires humidity/VPD and pressure contracts, explicit soil
moisture stress, new PFT parameters, and C3/C4 calibration. The first version
should use prescribed leaf temperature equal to air temperature; iterative
leaf-energy balance and plant hydraulics are later extensions.

## 6. Later extensions

- rotations, sequential and simultaneous crops, and shared soil resources;
- broader output/observation operators and global validation;
- gradient calibration, data assimilation, and hybrid ML processes;
- automatic spatial fallback batches, multi-GPU, and MPI;
- alternative soil hydraulic inputs and broader land-system coupling.
