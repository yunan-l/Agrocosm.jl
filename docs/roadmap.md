# Agrocosm.jl roadmap

This roadmap separates the current LPJmL-informed single-crop foundation from
the broader goal of a modular, differentiable, GPU-accelerated crop-model
framework. Completion means that a process is implemented, audited, and has an
appropriate CPU/GPU regression test; it does not by itself imply global
validation.

## Phase 1 — single-crop process foundation

### Completed or substantially completed

- C3/C4 photosynthesis, water-limited `lambda`, phenology, canopy growth,
  carbon allocation, respiration, sowing, harvest, and residue routing.
- Crop N allocation, demand, uptake, prescribed fertilizer/manure, and soil
  mineral-N supply.
- Five-layer soil water, snow, temperature, freeze--thaw, phase-change energy,
  and percolation enthalpy, using the current Phase-1 approximations.
- Non-methane soil C--N decomposition, mineralization, immobilization,
  nitrification, denitrification, volatilization, leaching, and fixed `c_shift`
  routing.
- LPJmL-informed full surface albedo and tillage--topsoil hydraulic coupling.
- Audited C3/C4 daily ordering for climate history, cultivation, albedo/PET,
  snow, soil preparation, C--N decomposition, crop processes, water removal,
  denitrification, and NH3 volatilization.
- Daily crop water-deficit output using LPJmL's `0--100%` definition.
- CPU/GPU process kernels, `Float32`/`Float64` support, mass/energy balance
  diagnostics, a high-level simulation API, and a 20-year rainfed-wheat
  notebook with numerically closed water, C, N, and energy ledgers.

### Phase 1 acceptance status

- The 20-year rainfed-wheat notebook has been rerun after the water-deficit
  output fix. Active-crop values are finite and physically bounded, and the
  water, C, N, and energy ledgers remain numerically closed.
- CPU/CUDA daily-output equivalence has been validated by the user.
- Public `save_checkpoint`/`restore_checkpoint!` APIs now write
  backend-independent files. An interrupted-file-restore trajectory is
  identical to an uninterrupted CPU simulation for crop, soil, output, and all
  four balance ledgers.
- The NH3 equation and units have been audited and accepted; further LPJmL
  input-parity comparison is not a Phase-1 requirement.

Phase 1 is complete. The rainfed-wheat notebook and short end-to-end tests
remain regression baselines rather than open implementation work.

## Phase 1 → Phase 2 architecture transition

The Terrarium-style separation of process configuration from numerical-state
lifecycle is now the Phase-2 architecture baseline:

- `ProcessModules` contains process choices and global parameters, with no
  evolving backend arrays.
- `ModelState` is the single canonical numerical tree, partitioned into
  `prognostic`, `fluxes`, `auxiliary`, `inputs`, `events`, `workspace`, and
  `output`, with crop and soil namespaces inside each lifecycle group.
- All crop, soil, climate, and balance process wrappers now select their arrays
  directly from `ModelState`; the bottom `@kernel` interfaces remain explicit
  array arguments.
- Checkpoint format v2 serializes prognostic state and restart-relevant inputs
  directly from the lifecycle tree. The temporary format-v1 compatibility
  path was removed together with the old runtime entry point.
- C3 and C4 old/new entry routes are exactly equal across every runtime array
  in the three-day regression (`1032/1032` assertions), and the complete CPU
  suite passes (`2533/2533`).

Next, define a one-day transition suitable for AD and select active parameters.

## Phase 2 — modular process alternatives and multi-crop architecture

The second phase turns the current process implementation into a framework for
scientific comparison and model development.

### Differentiable runtime foundation

1. Define the one-day transition boundary over `ProcessModules` and
   `ModelState`, keeping I/O and reporting outside the differentiated region.
2. Add an Enzyme CPU smoke test for selected active parameters and prognostic
   state, followed by a CUDA differentiation test where supported.
3. Audit sowing, harvest, fertilization, clamps, and iterative solvers and
   document which gradients are intentionally piecewise or inactive.
4. Add differentiability regressions without reintroducing domain-container
   aliases into the active runtime state.

### Spin-up, restart, and output completion

1. Implement soil/ecosystem spin-up for consistent initial C, N, water, and
   thermal states, including the post-spin-up `c_shift` routing configuration.
2. Validate restart continuity across the spin-up-to-transient boundary.
3. Complete the soil and climate time-series output chains and define stable
   output metadata without duplicating daily process calculations.
4. Build the independent `AgrocosmData.jl` input layer following the dedicated
   [data roadmap](agrocosm_data_roadmap.md): canonical `cellid` indexing,
   land-use/PFT masks, HWSD C/N preprocessing, and time-blocked global forcing.

### Interchangeable process modules

Define stable interfaces so a simulation can select among scientifically
documented alternatives without changing the main daily driver. The first
target is **photosynthesis**: retain the current LPJmL-informed C3/C4 pathway
as a reference implementation and add alternative C3/C4 photosynthesis or
stomatal-conductance formulations behind the same interface. The same pattern
can subsequently support alternative respiration, phenology, allocation,
soil-temperature, and decomposition schemes.

Each alternative must provide its assumptions, parameters, CPU/GPU behaviour,
numerical tests, and a comparison against the reference process. Process choice
must be explicit in simulation configuration and output metadata.

### Multi-crop stands and management

1. Flatten stand × crop into the existing batch dimension, keeping each
   crop--grid-cell instance independent and GPU-friendly.
2. Add stand/crop indexing and aggregation back to grid-cell outputs.
3. Introduce crop rotations, sequential crops, and later simultaneous crops
   with shared soil resources.
4. Extend management beyond prescribed events: irrigation, dynamic sowing, and
   more complete fertilizer/manure strategies.

## Phase 3 — differentiable, calibrated, and scalable simulations

- Replace or expose non-smooth numerical choices where required for robust
  automatic differentiation, while retaining scientifically meaningful
  reference modes.
- Build gradient-based parameter calibration and sensitivity-analysis workflows.
- Assimilate LAI, GPP, evapotranspiration, biomass, and yield observations.
- Develop hybrid process--machine-learning components with explicit physical
  constraints.
- Add reproducible CPU/GPU performance benchmarks for increasingly large
  spatial batches.

## Phase 4 — equilibrium, validation, and Earth-system coupling

- Establish multi-site and global validation protocols.
- Support high-resolution regional and large-domain GPU simulations.
- Couple Agrocosm, where useful, to broader land and Earth-system frameworks
  while retaining an independent crop-model API.