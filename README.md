# Agrocosm.jl

[![Build Status](https://github.com/yunan-l/Agrocosm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yunan-l/Agrocosm.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/yunan-l/Agrocosm.jl/actions/workflows/Documentation.yml/badge.svg)](https://yunan-l.github.io/Agrocosm.jl/dev/)

**🧑‍🌾 💧 ☀️ 🌾 🚀 Fast and flexible Julia framework for agricultural ecosystem modelling across scales.**

Agrocosm.jl is a framework for building a new generation of process-based crop models. It is a
crop–soil model with coupled water, carbon, nitrogen, and energy processes that runs on both CPUs and
GPUs, is easy to use and to extend, and is written in Julia to make physically based simulation,
differentiable programming, high-performance computing, and machine-learning workflows available
within one modelling environment.

Agrocosm is a **downstream package built on the [Terrarium.jl](https://github.com/NumericalEarth/Terrarium.jl)
land-modelling framework**. Terrarium supplies all of the infrastructure — grids, state, continuous-time
timestepping, CPU/GPU architectures, parameters, I/O, checkpointing, and differentiability — together
with the physical soil and surface processes (soil energy, Richards hydrology, the surface energy
balance). Agrocosm contributes the crop-specific parts: C3/C4 photosynthesis, carbon and nitrogen
allocation, phenology, a coupled soil carbon–nitrogen biogeochemistry, crop management, and a
top-level managed-crop land model.

The crop physiology takes the crop module of [LPJmL](https://github.com/PIK-LPJmL/LPJmL) as a
scientific reference. Agrocosm is **not** a line-by-line port of LPJmL: it re-expresses the relevant
process logic as continuous-time, GPU-ready, differentiable Terrarium processes.

Read the [online documentation](https://yunan-l.github.io/Agrocosm.jl/dev/) for model concepts, the
process API, CPU/GPU execution, differentiability, and the API reference.

## Vision

We want Agrocosm to be:

- **Fully GPU-compatible**, from a single site to large ensembles of grid cells
- **Differentiable**, enabling gradient-based calibration, sensitivity analysis, data assimilation,
  and hybrid modelling
- **Process-based and auditable**, with explicit carbon, nitrogen, water, and energy dynamics
- **Modular and extensible**, comparing alternative process representations without rebuilding the
  full model
- **Open and community-oriented**, providing a foundation that can incorporate new crop physiology and
  collaborate with the wider crop-modelling community.

## Architecture

Agrocosm assembles a Terrarium `LandModel` from Terrarium's physical soil and surface with Agrocosm's
crop processes. The division of labour:

| Layer | Provided by | Contents |
| --- | --- | --- |
| Infrastructure | **Terrarium** | grids, state, continuous-time timesteppers, CPU/GPU/Reactant architectures, parameters, I/O, checkpointing, differentiability |
| Physical soil & surface | **Terrarium** | soil energy, Richards hydrology, stratigraphy, surface energy balance, atmosphere coupling |
| Crop physiology | **Agrocosm** | C3/C4 photosynthesis + λ solver, stomatal conductance, phenology (prognostic heat units), carbon & nitrogen pools with the N→Vcmax feedback, root distribution, plant-available water |
| Soil biogeochemistry | **Agrocosm** | `CropSoilBiogeochemistry`: prognostic litter/fast/slow carbon and ammonium/nitrate, with mineralization, nitrification, denitrification, and the plant↔soil flux loop |
| Management | **Agrocosm** | sowing/harvest as continuous-time-sanctioned discrete callbacks; fertilizer as a continuous input flux |

All crop processes are continuous-time `AbstractProcess`/`AbstractVegetation`/`AbstractSoilBiogeochemistry`
implementations written as type-stable, allocation-free, throw-free kernels, so they run on CPU and GPU
and differentiate with [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) (directly on the CPU and
through [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl)).

## Installation

Agrocosm is not yet registered in the Julia General registry. Clone the source and instantiate its
project environment:

```bash
git clone https://github.com/yunan-l/Agrocosm.jl.git
cd Agrocosm.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Agrocosm tracks Terrarium via a `[sources]` entry in `Project.toml`, so `Pkg.instantiate()` fetches
the required Terrarium revision automatically. An NVIDIA GPU with a working
[CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) install is needed only for GPU execution and the
Reactant tests.

## Quick start

Build a managed-crop model for a crop functional type and run it on a soil column:

```julia
using Agrocosm
using Terrarium

# A single soil column (20 layers, exponentially stretched).
grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# A maize (C4) crop model: crop vegetation + crop soil carbon–nitrogen biogeochemistry.
model = CropModel(grid, crop_pft("maize"))

integrator = initialize(model; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)
set!(integrator.state.surface_shortwave_down, 400.0)

run!(integrator; steps = 100, Δt = 600.0)

interior(integrator.state.leaf_area_index)          # crop canopy
interior(integrator.state.gross_primary_production)  # crop GPP
interior(integrator.state.soil_nitrate)              # soil mineral nitrogen
```

The 12 LPJmL crop functional types are available through `crop_pft(name_or_number)` (e.g.
`crop_pft("temperate cereals")`, `crop_pft(3)`) or the `cft1`…`cft12` presets.

## Crop management

Sowing and harvest are genuine discrete lifecycle events, implemented as
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) callbacks on a `Simulation` (the
continuous-time framework's sanctioned exception); fertilizer is a continuous input flux to the soil
mineral nitrogen:

```julia
using Oceananigans: TimeInterval

simulation = Simulation(integrator; Δt = 600.0, stop_time = 200 * 24 * 3600.0)

calendar = CropCalendar(Float64; sowing_day = 120, harvest_day = 280, residue_fraction = 0.25)
add_crop_management!(simulation, calendar)   # seeds the stand at sowing; exports grain + returns residue at harvest

fertilization = CropFertilization(Float64; application_rate = 1e-7, application_start_day = 120, application_end_day = 160)
add_crop_fertilization!(simulation, fertilization)   # continuous fertilizer over the window

run!(simulation)
```

## Differentiability

Because the crop processes are differentiable Terrarium kernels, you can take derivatives straight
through a model integration — the basis for gradient-based calibration and hybrid physics/ML models.
Reverse-mode adjoints of the crop scalar primitives are checked on the CPU in the test suite
(`test/crop/test_differentiability.jl`); worked model-level examples live in `examples/autodiff/`:

- `crop_soil_reactant.jl` — compiling and running the model through Reactant.
- `differentiating_crop_soil_reactant.jl` — reverse-mode AD of the compiled rollout (Enzyme + Reactant).

## Testing

Run the CPU test suite (includes the crop physics and CPU Enzyme adjoint checks):

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The CPU-vs-Reactant correctness tests and the Reactant reverse-mode autodiff tests live in a dedicated
environment and require Julia 1.12 with a working Reactant/CUDA install:

```bash
julia +1.12 --project=test/reactant test/reactant/runtests.jl
```

## Contributing

Contributions, ideas, issue reports, and process-comparison experiments are very welcome. Agrocosm is
most useful when crop physiologists, Earth system modellers, numerical scientists, and
machine-learning researchers can inspect, question, and improve its assumptions together. Framework
conventions live in [`AGENTS.md`](AGENTS.md) and Terrarium's `AGENTS.md`; please open a GitHub issue to
start a discussion.

## Acknowledgements

Agrocosm.jl is a research project developed with the support of the
[Earth System Modeling group](https://www.asg.ed.tum.de/esm/home/) at the Technical University of
Munich (TUM) and the
[FutureLab on Artificial Intelligence](https://www.pik-potsdam.de/en/institute/departments/complexity-science/research/artificial-intelligence)
at the Potsdam Institute for Climate Impact Research (PIK). The author acknowledges funding from the
China Scholarship Council (grant agreement no. 202303250017) and the Horizon Europe ClimTip project
(grant agreement no. 101137601).

## License

Agrocosm.jl is released under the [European Union Public Licence v1.2](https://eupl.eu/1.2/en). You are
encouraged to copy, modify, and build upon our code to advance your research.
