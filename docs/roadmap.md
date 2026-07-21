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
  and percolation enthalpy.
- Non-methane soil C--N decomposition, mineralization, immobilization,
  nitrification, denitrification, volatilization, leaching, and fixed `c_shift`
  routing.
- CPU/GPU process kernels, `Float32`/`Float64` support, mass/energy balance
  diagnostics, and a high-level simulation API.

### Current closing work

1. Confirm the minimal checkpoint boundary and the shared post-spin-up
   `c_shift` configuration on CPU and GPU.
2. Maintain a short end-to-end rainfed-wheat reference simulation as a
   regression baseline for the single-crop model.
3. Audit every daily process and its ordering against the relevant LPJmL source
   path; record intentional differences rather than pursuing impractical
   one-to-one output matching.
4. Improve public examples, input documentation, and benchmark reporting.

## Phase 2 — modular process alternatives and multi-crop architecture

The second phase turns the current process implementation into a framework for
scientific comparison and model development.

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

- Implement soil/crop spin-up and restart workflows for consistent equilibrium
  initial conditions.
- Establish multi-site and global validation protocols.
- Support high-resolution regional and large-domain GPU simulations.
- Couple Agrocosm, where useful, to broader land and Earth-system frameworks
  while retaining an independent crop-model API.

## Deliberately deferred scope

Methane/wetland pathways and waterlogging stress are not part of the current
annual-crop core. They should be considered only when a rice/wetland or
water-saturated crop use case requires them.
