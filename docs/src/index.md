# Agrocosm.jl

**Agrocosm.jl is a process-based, GPU-capable, differentiable Julia model of crop–soil water, carbon,
nitrogen, and energy dynamics.**

Agrocosm is a **downstream package built on the [Terrarium.jl](https://github.com/NumericalEarth/Terrarium.jl)
land-modelling framework**. Terrarium supplies the infrastructure — grids, state, continuous-time
timestepping, CPU/GPU/Reactant architectures, parameters, I/O, checkpointing, and differentiability —
and the physical soil and surface processes (soil energy, Richards hydrology, the surface energy
balance). Agrocosm contributes the crop-specific parts:

- **Crop physiology** — C3/C4 photosynthesis with the λ water-coupling solver, stomatal conductance,
  phenology driven by prognostic heat units, carbon and nitrogen pools with the nitrogen→Vcmax
  feedback, root distribution, and plant-available water.
- **Soil biogeochemistry** — `CropSoilBiogeochemistry`, with prognostic litter/fast/slow soil carbon
  and ammonium/nitrate mineral nitrogen, coupled by mineralization, nitrification, and denitrification,
  and closing the plant↔soil flux loop (crop nitrogen uptake and litter return).
- **Management** — sowing and harvest as discrete Oceananigans callbacks (the continuous-time
  framework's sanctioned exception) and fertilizer as a continuous input flux.
- **A managed-crop model** — `CropModel`, assembling all of the above into a Terrarium `LandModel`, and
  the 12 LPJmL crop functional types.

The crop physiology uses [LPJmL](https://github.com/PIK-LPJmL/LPJmL) as a scientific reference, but is
an independent Julia implementation re-expressed as continuous-time, differentiable Terrarium
processes rather than a line-by-line translation.

## Why Agrocosm?

- One process implementation runs on CPU and GPU (and compiles through Reactant/XLA).
- Both `Float32` and `Float64` are supported.
- The crop processes differentiate with Enzyme — directly on the CPU and through Reactant — for
  gradient-based calibration and hybrid physics/ML models.
- Crop physics is modular: alternative photosynthesis, phenology, or biogeochemistry can be swapped in
  without rebuilding the model.

Start with [Getting started](@ref) to build and run a crop model, then see the [API reference](@ref).

```@contents
Pages = ["getting_started.md", "api.md"]
Depth = 2
```
