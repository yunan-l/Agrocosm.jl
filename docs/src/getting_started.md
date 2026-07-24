# Getting started

## Installation

Agrocosm is not yet registered in Julia's General registry. Clone the repository and instantiate its
project environment:

```bash
git clone https://github.com/yunan-l/Agrocosm.jl.git
cd Agrocosm.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Agrocosm tracks [Terrarium.jl](https://github.com/NumericalEarth/Terrarium.jl) via a `[sources]` entry
in `Project.toml`, so instantiation fetches the required Terrarium revision automatically. A CUDA
device is optional and needed only for GPU execution and the Reactant tests (see below).

## Build and run a crop model

Agrocosm's entry point is [`CropModel`](@ref), which assembles a Terrarium `LandModel` from the crop
vegetation and the crop soil carbon–nitrogen biogeochemistry for a chosen crop functional type:

```julia
using Agrocosm
using Terrarium

# A single soil column, 20 exponentially stretched layers.
grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# Maize (CFT 3, a C4 crop).
model = CropModel(grid, crop_pft("maize"))

integrator = initialize(model; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)
set!(integrator.state.surface_shortwave_down, 400.0)

run!(integrator; steps = 100, Δt = 600.0)
```

The 12 LPJmL crop functional types are addressed by name or number through `crop_pft`, or via the
`cft1`…`cft12` presets:

```julia
crop_pft("temperate cereals")   # CFT 1
crop_pft(3)                      # maize
```

## Inspect results

State variables are Terrarium `Field`s; read their interiors with `interior`:

```julia
lai = interior(integrator.state.leaf_area_index)              # crop canopy
gpp = interior(integrator.state.gross_primary_production)     # crop GPP
biomass = interior(integrator.state.crop_biomass)             # crop carbon pool
nitrate = interior(integrator.state.soil_nitrate)             # soil mineral nitrogen
litter = interior(integrator.state.litter_carbon)             # soil litter carbon
```

## Crop management

Sowing and harvest are discrete lifecycle events, applied through Oceananigans callbacks on a
`Simulation`; fertilizer is a continuous input flux. See [`CropCalendar`](@ref),
[`add_crop_management!`](@ref), and [`add_crop_fertilization!`](@ref):

```julia
simulation = Simulation(integrator; Δt = 600.0, stop_time = 200 * 24 * 3600.0)

calendar = CropCalendar(Float64; sowing_day = 120, harvest_day = 280, residue_fraction = 0.25)
add_crop_management!(simulation, calendar)

fertilization = CropFertilization(Float64; application_rate = 1e-7, application_start_day = 120, application_end_day = 160)
add_crop_fertilization!(simulation, fertilization)

run!(simulation)
```

At `sowing_day` the stand is seeded and the phenological clock is reset; at `harvest_day` the grain is
exported, the residue is returned to the soil litter over the root zone (mass-conserving), and the
stand is cleared.

## Differentiability

The crop processes differentiate with Enzyme, on the CPU and through Reactant. Reverse-mode adjoints of
the crop scalar primitives are checked on the CPU in `test/crop/test_differentiability.jl`; model-level
examples are in `examples/autodiff/`:

- `crop_soil_reactant.jl` — compiling and running the model through Reactant.
- `differentiating_crop_soil_reactant.jl` — reverse-mode AD of the Reactant-compiled rollout.

The CPU-vs-Reactant correctness tests and the Reactant autodiff tests live in `test/reactant/` and run
on Julia 1.12:

```bash
julia +1.12 --project=test/reactant test/reactant/runtests.jl
```
