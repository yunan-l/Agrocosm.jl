# Getting started

## Installation

Agrocosm is not yet registered in Julia's General registry. Install it from
GitHub or clone the repository:

```julia
import Pkg
Pkg.add(url = "https://github.com/yunan-l/Agrocosm.jl")
```

For development:

```bash
git clone https://github.com/yunan-l/Agrocosm.jl.git
cd Agrocosm.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Julia 1.10 is currently supported. A CUDA device is optional.

## Run the included wheat example

The repository includes initial conditions and ten years of daily forcing in
`examples/`. From the repository root:

```julia
using Agrocosm
using JLD2

initial_data = load("examples/initial_wheat.jld2", "initial_data")
climate = load("examples/climate_2000_2009.jld2", "climate")

simulation = initialize_simulation(
    cft1,
    initial_data;
    indices = [1],
    device = identity,
    T = Float32,
    days = size(climate.temp, 1),
    auto_fertilizer = false,
)

run_simulation!(simulation, climate)
summary = simulation_summary(simulation)
```

`cft1` is the current wheat parameter set. Use `cft3` for the available C4
pathway. Parameter sets are research defaults, not universal cultivar
calibrations.

## Inspect results

Completed daily outputs are stored under `simulation.output`:

```julia
npp = Array(simulation.output.crop.npp)
lai = Array(simulation.output.crop.lai)
soil_water = Array(simulation.output.soil.water_storage)
water_deficit_percent = Array(simulation.output.crop.water_deficit)
```

The canonical numerical state is under `simulation.state`, for example:

```julia
leaf_carbon = simulation.state.prognostic.crop.carbon.leaf
soil_liquid_water = simulation.state.prognostic.soil.water.storage
daily_percolation = simulation.state.fluxes.soil.water.percolation
```

Do not treat daily fluxes and auxiliary fields as restart state. See
[State variables](@ref).
