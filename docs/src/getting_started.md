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

## Driving the model with climate data

Instead of constant forcing, drive the model with tabulated climate series through Terrarium's
time-varying input sources. [`surface_climate_inputs`](@ref) packs each daily series into an
Oceananigans `FieldTimeSeries` wrapped in a `FieldTimeSeriesInputSource`; passing the result as
`inputs` lets `run!` interpolate the forcing to the current time on every step — no manual per-step
`set!`:

```julia
using JLD2
climate = load("examples/climate_2000_2009.jld2", "climate")   # 3650 days × cells
temp = Float64.(climate.temp[:, 1])                            # °C
shortwave = Float64.(climate.swdown[:, 1])                     # W/m²
times = (0:length(temp) - 1) .* Terrarium.seconds_per_day(Float64)

inputs = surface_climate_inputs(grid, times; air_temperature = temp, surface_shortwave_down = shortwave)
integrator = initialize(model; inputs, initializers = (temperature = 2.0,))
run!(integrator; steps = 144, Δt = 600.0)                     # one day of forcing
```

`surface_climate_inputs` also accepts a `time × column` matrix per variable, giving each column of a
multi-column / global [`ColumnRingGrid`](@extref Terrarium.ColumnRingGrid) its own series. The complete
worked runs are in `examples/`:

- `wheat_gpp_npp.jl` — a ten-year single-column wheat run driven by the daily climate.
- `wheat_gpp_npp_global.jl` — the same over a batch of columns on a global grid.

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
