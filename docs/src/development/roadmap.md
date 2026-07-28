# Roadmap

## Completed foundation

Phase 1 is complete. Agrocosm has an LPJmL-informed single-crop C3/C4 process
baseline, CPU/GPU kernels, lifecycle-organized numerical state, daily balance
diagnostics, checkpoints, streamed output, a high-level simulation API, and a
public one-day transition. The scientific documentation is organized as a
Model processes overview with dedicated crop, soil, climate/surface,
numerics, and initialization/warm-up pages.

Launch hardening is also complete at the CPU regression level. Checkpoints now
bind compact cell identity and PFT identity; production warm-up has a strict
convergence gate and explicit memory accounting; process, initialization, and
output updates use backend kernels with synchronization at lifecycle
boundaries rather than after every kernel. Legacy `_reference!` paths have
been removed. The current CPU suite passes 2092 tests.

The AgrocosmData core is also substantially complete:

- canonical `cellid` indexing on the `720 × 280` grid;
- the 12-crop registry and explicit management-band maps;
- compact soil and management readers;
- HWSD 2.x SOC/total-N preprocessing and native initial-state construction;
- bounded climate-block reading, calendar/unit normalization, annual CO₂, and
  one-block prefetch;
- model-facing `model_initial_data` and `climate_forcings` adapters.
- a configuration-driven server utility that extracts one rainfed wheat band
  from 2015 multi-PFT management data and the 2015–2016 daily forcing.

This means new grids no longer require an LPJmL restart. It does not yet mean
that the global production workflow is complete.

## Current phase: global production readiness

Work in this phase is ordered as follows:

1. Use 2015 rainfed-wheat `landfrac > 0` only to select the fixed compact cell
   set. Land fraction is not a multiplier in crop or soil process equations;
   all other management inputs are likewise fixed at their 2015 values.
2. The local ten-cell HWSD + 2015–2016 forcing smoke test and restartable
   streamed `agricultural_warmup!` are complete. Native post-warm-up and 2015
   boundary checkpoint/restart are now part of the production runner.
3. The canonical-grid HWSD pipeline and QC/fallback contracts are implemented;
   retain the server product and its QC report with every production run.
4. Submit the current bounded-memory CPU and single-GPU workflows through
   Slurm rather than interactive nodes. Both jobs run their complete backend
   regression suite before entering production. CPU and GPU outputs and
   checkpoints must use separate directories.
5. Validate memory, throughput, grid reconstruction, finite/non-negative state,
   CPU/GPU agreement, sampled C/N/water/energy closure, and 2015→2016 native
   checkpoint continuity. The latest unified GPU suite still requires a fresh
   server run after its isolated-test entry point was repaired.
6. Use target-constrained warm-up with a strict production gate. A previous
   100-year global diagnostic reached about 96.15% converged cells, so the
   current code must identify and review the remaining cells rather than
   silently writing a production checkpoint.
7. Retain the interim HWSD 40:60 fast/slow split as the reproducible first-run
   baseline, but keep it under review. The real ten-cell warm-up remains
   transient after ten years; changing only the initial ratio would not resolve
   the continuing total C/N decline.

The phase is complete when the current commit passes both full backend suites
and a 2015-selected global crop domain can be initialized from native data,
warmed under the declared convergence contract, checkpointed, run across a
year boundary, and reconstructed to the canonical grid without LPJmL-derived
state.

## Phase 2: differentiable transition

1. Declare the active parameter/state boundary for the existing one-day
   transition.
2. Add Enzyme CPU smoke tests, finite-difference gradient checks, and explicit
   policies for discrete sowing, harvest, fertilization, clamps, and failed
   crops.
3. Add CUDA differentiation only after the CPU gradient path is stable.
4. Keep data loading, warm-up, checkpoints, and reporting outside the
   differentiated region.

## Phase 3: alternative processes and multi-crop operation

- Retain the current LPJmL-informed canopy exchange as the reference pathway.
- Add an optional Farquhar photosynthesis + Medlyn conductance + simplified
  Penman–Monteith pathway after humidity/VPD and pressure forcing contracts are
  available. Compare it against the reference before changing defaults.
- Add crop rotations, sequential crops, stand/crop indexing, and shared-soil
  management.
- Complete broader soil/climate outputs as required by validation workflows.

## Later work

Later phases cover gradient-based calibration, data assimilation, hybrid
process–machine-learning components, multi-site/global validation, spatial
fallback batching, multi-GPU/MPI execution, and coupling to broader land or
Earth-system frameworks.
