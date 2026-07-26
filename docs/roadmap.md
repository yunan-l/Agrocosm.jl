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
  ledgers untouched while retaining warmed state.

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
- configuration-driven extraction of a single rainfed-wheat management band
  and the first two 365-day climate years for a bounded global test dataset.

Remaining data-layer work is production hardening rather than new loader
architecture:

- execute and quality-control the full canonical-grid HWSD product;
- preserve full source/provenance manifests for server runs;
- benchmark real server NetCDF access and add a canonical cache only if direct
  compact reads are too slow.

Annual crop activation, warm-up, backend transfer, state evolution, and global
execution remain responsibilities of Agrocosm.jl, not AgrocosmData.jl.

## 3. Immediate production sequence

### 3.1 Annual land-use activation

- Pass the fixed-union `CropMask.selection` into initialization.
- Add annual `active` and crop-fraction inputs to runtime state.
- Gate cultivation, fertilizer/manure, crop uptake, and crop output by annual
  activity while continuing soil water, heat, and C/N processes in fallow
  cells.
- Test zero-to-positive, positive-to-zero, and continuously active sequences.

Acceptance: fixed allocation and dynamically active execution reproduce
separate active-year reference runs without reallocating backend state.

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

### 3.3 Global rainfed-wheat smoke test

- Generate the two-year rainfed-wheat subset on the server and retain a source
  manifest; generate the complete canonical-grid HWSD product and review its
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
