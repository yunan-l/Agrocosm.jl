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
- Checkpoint schema validation for compact `cell_ids` and CFT identity, strict
  warm-up convergence gating, and explicit warm-up/cached-forcing memory
  accounting.
- Backend kernels throughout daily processes and initialization/output state
  copies, with synchronization at initialization and daily lifecycle
  boundaries instead of after each kernel. Legacy `_reference!` paths are
  removed.

This foundation remains the scientific regression baseline. Alternative
processes must demonstrate their differences against it rather than silently
replacing it.

## 2. AgrocosmData status

Milestones 1–5 are substantially complete at the code and fixture-test level:

- dataset catalog and versioned backend-neutral contracts;
- canonical grid selection and compact/global round trips;
- 12-CFT registry, explicit 64/32/24/16-band mappings, and crop masks;
- soil-code properties, pH, sowing date, PHU, fertilizer, manure, residue, and
  land-use readers;
- HWSD SOC/total-N aggregation, vertical remapping, uncertainty/fallback
  provenance, field-capacity water, and native initial state;
- daily temperature, precipitation, net longwave, and downward shortwave
  streaming; annual CO₂ alignment; 365-day normalization; block prefetch;
- full ten-cell equivalence through `model_initial_data` and
  `climate_forcings`.
- configuration-driven extraction of selected CFT management fields and
  climate periods without materializing complete source datasets.

Remaining data-layer work is production hardening rather than new loader
architecture:

- validate canonical-grid HWSD coverage and stock conservation;
- preserve source/provenance manifests for generated data products;
- benchmark representative NetCDF access and add a canonical cache only if direct
  compact reads are too slow.

Warm-up, backend transfer, state evolution, and global
execution remain responsibilities of Agrocosm.jl, not AgrocosmData.jl.

## 3. Simulation acceptance criteria

### 3.1 Configured crop domain

- Select compact cells from the configured CFT, water system, years and land-use mask.
- Use `landfrac` only for selection and provenance; do not multiply any model
  process or reported crop quantity by fractional area.
- Resolve sowing date, PHU, fertilizer, manure and residue settings through
  the configured fixed or transient management policy.

Acceptance: every selected cell runs one crop stand in stable compact
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
then checks checkpoint/restart at a production-year boundary.

### 3.3 Backend and scale validation

- Validate initialization and cross-year state propagation on a bounded fixture.
- Review the canonical-grid HWSD product's
  coverage and conservation summaries.
- Compare CPU and GPU execution using identical compact cell ordering.
- Stream climate and monthly/annual output; avoid full daily global ledgers.
- Check NaN/Inf, invalid negative pools, crop lifecycle failures, memory peak,
  throughput, restart continuity, and sampled or online balance closure.
- Reconstruct outputs to `720 × 280` by `cellid` and verify mask alignment.
- Continue through the second forcing year to test cross-year state and
  checkpoint/restart continuity after the first-year smoke test passes.

Acceptance: the complete year finishes within estimated memory, CPU/GPU
differences meet declared tolerances, and the second-year restart/reassembly is
deterministic.

Experiment-specific launch plans, convergence summaries and acceptance records
are maintained outside the package. A completed calibration alone does not
establish suitability for every production configuration.

### 3.4 HWSD pool-allocation decision

The current native initialization conserves HWSD layer SOC and total N using a
documented 40:60 fast/slow split and zero litter. Do not add an elaborate
equilibrium allocator without evidence. Inspect configured warm-up diagnostics for
initial respiration pulses, litter/fast-pool stabilization, mineral-N drift,
and total C/N trajectories. If needed, implement a constrained allocation that
preserves every layer total and records uncertainty.

The 40:60 initialization is not an equilibrium claim. Review total-stock drift
as well as pool fractions; changing only the initial ratio does not establish
equilibrium. Allocation and convergence decisions belong to the declared
initialization contract of each experiment.

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
moisture stress, new CFT parameters, and C3/C4 calibration. The first version
should use prescribed leaf temperature equal to air temperature; iterative
leaf-energy balance and plant hydraulics are later extensions.

## 6. CFT patch production sequence

- Keep CFT as the public model/data contract. Legacy LPJmL `pft` is permitted
  only while reading source dimensions and historical allocation files.
- Each `(cft_id, irrigated)` patch has independent crop, soil, litter, and
  management state. `landfrac` selects patches and weights reconstructed
  outputs; it does not scale cell-level process equations.
- Derive and store one soil-pool allocation product per selected patch batch.
- Verify each batch on CPU/GPU, native checkpoint/restart, sampled balance
  closure, and canonical-grid aggregation before adding rotations.

## 7. Later extensions

- rotations, sequential and simultaneous crops, and shared soil resources;
- broader output/observation operators and global validation;
- gradient calibration, data assimilation, and hybrid ML processes;
- automatic spatial fallback batches, multi-GPU, and MPI;
- alternative soil hydraulic inputs and broader land-system coupling.
