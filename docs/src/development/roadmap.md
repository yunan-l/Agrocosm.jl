# Roadmap

## Completed foundation

Phase 1 is complete. Agrocosm has an LPJmL-informed single-crop C3/C4 process
baseline, CPU/GPU kernels, lifecycle-organized numerical state, daily balance
diagnostics, checkpoints, streamed output, a high-level simulation API, and a
public one-day transition. The scientific documentation is organized as a
Model processes overview with dedicated crop, soil, climate/surface,
numerics, and initialization/warm-up pages.

Launch hardening is also complete at the CPU regression level. Checkpoints now
bind compact cell identity and CFT identity; production warm-up has a strict
convergence gate and explicit memory accounting; process, initialization, and
output updates use backend kernels with synchronization at lifecycle
boundaries rather than after every kernel. Legacy `_reference!` paths have
been removed. Process and diagnostics tests now exercise the same `ModelState`
interface as the runtime. Validation scope is documented separately from
experiment-specific test results.

The AgrocosmData core is also substantially complete:

- canonical `cellid` indexing on the `720 × 280` grid;
- the 12-crop registry and explicit management-band maps;
- compact soil and management readers;
- HWSD 2.x SOC/total-N preprocessing and native initial-state construction;
- bounded climate-block reading, calendar/unit normalization, annual CO₂, and
  one-block prefetch;
- model-facing `model_initial_data` and `climate_forcings` adapters.
- configuration-driven extraction of selected CFT management bands and
  climate periods using bounded input reads.

This means new grids no longer require an LPJmL restart. It does not yet mean
that the global production workflow is complete.

## Production validation requirements

- Derive compact cells from the configured crop mask. Land fraction does not
  multiply crop or soil process equations.
- Preserve data provenance, uncertainty and quality-control reports for
  initialization products.
- Verify bounded memory, grid reconstruction, finite/non-negative state,
  CPU/GPU agreement, balance closure and cross-year checkpoint continuity.
- Apply the declared warm-up allocation and convergence criteria; a completed
  warm-up is not automatically an equilibrium initial state.
- Validate each selected CFT and water-system configuration on the relevant
  backend before relying on its output scientifically.

Experiment periods, management scenarios, scheduler settings, run progress and
case-specific results are maintained outside the public package documentation.

## Differentiable transition

The optional Enzyme extension supports parameter objectives and fixed-event
C3/C4 weather-to-yield sensitivities on CPU. Development priorities include
broader gradient regression coverage and explicit treatment of nonsmooth
events. Keep data loading, warm-up, checkpoints and reporting outside the
differentiated region. Full-model CUDA differentiation requires separate
validation; CPU/GPU forcing-copy equivalence alone does not establish it.

## Later production extensions

- The public interface is CFT-based (`CFTParameters`, `CFTRegistry`, `cft_id`)
  while legacy LPJmL `pft` identifiers are accepted only while decoding source
  NetCDF or pre-migration allocation files.
- Complete the independent-soil 12-CFT × rainfed/irrigated patch batches,
  per-batch soil-pool allocation products, and landfrac-weighted output
  aggregation. CFT 1 rainfed is the reference regression case.
- Retain the current LPJmL-informed canopy exchange as the reference pathway.
- Add an optional Farquhar photosynthesis + Medlyn conductance + simplified
  Penman–Monteith pathway after humidity/VPD and pressure forcing contracts are
  available. Compare it against the reference before changing defaults.
- Add crop rotations and sequential crops after the independent-patch workflow
  is validated; do not introduce shared soil state by default.
- Complete broader soil/climate outputs as required by validation workflows.

## Later work

Later phases cover gradient-based calibration, data assimilation, hybrid
process–machine-learning components, multi-site/global validation, spatial
fallback batching, multi-GPU/MPI execution, and coupling to broader land or
Earth-system frameworks.
