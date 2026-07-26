# Roadmap

## Completed foundation

Phase 1 is complete. Agrocosm has an LPJmL-informed single-crop C3/C4 process
baseline, CPU/GPU kernels, lifecycle-organized numerical state, daily balance
diagnostics, checkpoints, streamed output, a high-level simulation API, and a
public one-day transition. The scientific documentation is organized as a
Model processes overview with dedicated crop, soil, climate/surface,
numerics, and initialization/warm-up pages.

The AgrocosmData core is also substantially complete:

- canonical `cellid` indexing on the `720 × 280` grid;
- the 12-crop registry and explicit management-band maps;
- compact soil and management readers;
- HWSD 2.x SOC/total-N preprocessing and native initial-state construction;
- bounded climate-block reading, calendar/unit normalization, annual CO₂, and
  one-block prefetch;
- model-facing `model_initial_data` and `climate_forcings` adapters.
- a configuration-driven server utility that extracts one rainfed wheat band
  from multi-PFT management data and the first two 365-day forcing years.

This means new grids no longer require an LPJmL restart. It does not yet mean
that the global production workflow is complete.

## Current phase: global production readiness

Work in this phase is ordered as follows:

1. Connect `CropMask.active` and annual crop fraction to the daily runtime.
   The allocation mask remains fixed while soil state continues through fallow
   years and crop cultivation is disabled where land use is inactive.
2. Extend `agricultural_warmup!` to consume restartable streamed climate
   blocks, retain its production-output isolation, and write a native warm-up
   checkpoint.
3. Generate the two-year global rainfed-wheat subset and the complete
   canonical-grid HWSD product on the server, retaining source manifests and
   HWSD coverage/conservation quality-control reports.
4. Run a one-year global smoke test over every land-use-selected cell: CPU
   first, then one GPU, using bounded forcing and output blocks. Use the second
   year to verify cross-year state and checkpoint/restart continuity.
5. Validate memory, throughput, grid reconstruction, finite/non-negative state,
   CPU/GPU agreement, and sampled or online C/N/water/energy closure.
6. Use the warm-up drift report to evaluate the interim HWSD 40:60 fast/slow
   split. Add constrained pool allocation only if the evidence shows that the
   fixed split creates unacceptable transients.

The phase is complete when a global crop mask can be initialized from native
data, warmed, checkpointed, run across a year boundary, and reconstructed to
the canonical grid without LPJmL-derived state.

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
