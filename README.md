# Agrocosm.jl

[![Build Status](https://github.com/yunan-l/Agrocosm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yunan-l/Agrocosm.jl/actions/workflows/CI.yml?query=branch%3Amain)

**A process-based crop modelling framework for differentiable, GPU-accelerated agricultural simulation in Julia.**

Agrocosm.jl is a framework for building a new generation of
process-based crop models. Technically it is a crop-soil model with water, carbon, nitrogen, and energy processes with a numerical design that can run on both CPUs and GPUs. It is easy to use and easy to extend. It is written in Julia to make physically based
simulation, differentiable programming, high-performance computing, and
machine-learning workflows available within one modelling environment.

The present implementation takes the crop module of
[LPJmL](https://github.com/PIK-LPJmL/LPJmL) as an scientific
reference. Agrocosm is **not** a line-by-line port of LPJmL. It is an
independent Julia implementation that preserves relevant process
logic while developing a GPU-aware and increasingly
differentiable model architecture.

> [!WARNING]
> Agrocosm.jl is under active development, but almost done as a standalone model.

## Vision

We want Agrocosm to be:

- **Fully GPU-compatible**, from a single site to large ensembles of grid cells
- **Differentiable**, enabling gradient-based calibration,
  sensitivity analysis, data assimilation, and hybrid modelling
- **Process-based and auditable**, with explicit carbon, nitrogen, water, and
  energy balance diagnostics
- **Modular and extensible**, comparing alternative process representations
  without rebuilding the full model
- **Open and community-oriented**, providing a foundation that can incorporate
  new crop physiology and collaborate with the wider crop-modelling community.

Agrocosm keeps an independent crop-model architecture, while aiming to remain suitable for
future coupling with land and Earth system modelling frameworks.

## Current scope

Agrocosm currently focuses on daily, gridded simulations of a single crop.

<table>
  <thead>
    <tr>
      <th>Model component</th>
      <th>Subsystem</th>
      <th>Current implementation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="2" align="center" valign="middle"><strong>Crop</strong></td>
      <td>Carbon and phenology</td>
      <td>C3/C4 photosynthesis, canopy growth, phenology, carbon allocation, and respiration.</td>
    </tr>
    <tr>
      <td>Nitrogen and management</td>
      <td>Crop N allocation and uptake; cultivation, fertilizer, harvest, and residue transfer.</td>
    </tr>
    <tr>
      <td rowspan="2" align="center" valign="middle"><strong>Soil</strong></td>
      <td>Physics</td>
      <td>Five-layer soil water, snow, temperature, freeze--thaw, and water/energy transport.</td>
    </tr>
    <tr>
      <td>Biogeochemistry</td>
      <td>Litter routing; coupled soil C--N decomposition; mineralization, immobilization, nitrification, denitrification, volatilization, and leaching.</td>
    </tr>
    <tr>
      <td align="center" valign="middle"><strong>Numerics and verification</strong></td>
      <td>Backend and precision</td>
      <td>CPU/GPU kernels via <a href="https://github.com/JuliaGPU/KernelAbstractions.jl">KernelAbstractions.jl</a>, <code>Float32</code>/<code>Float64</code> support.</td>
    </tr>
  </tbody>
</table>

## What Agrocosm is not yet

The following are planned, but not yet part of the current model:

- **Near-term:** ecosystem spin-up, a complete nitrogen--photosynthesis
  feedback, and broader management processes.
- **Next model generation:** multi-crop stands, rotations, and dynamic sowing.
- **Long-term:** end-to-end differentiable process pathways, global validation,
  and hybrid process--machine-learning applications.

Methane and waterlogging stress are intentionally outside the current scope.

## Installation

Agrocosm is not yet registered in the Julia General registry. Clone the source
and instantiate its project environment:

```bash
git clone https://github.com/yunan-l/Agrocosm.jl.git
cd Agrocosm.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

A NVIDIA GPU and a working [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) installation are needed only for GPU execution and GPU tests.

## Quick start

```julia
using Agrocosm

# `initial_data` contains soil initial states.
# `climate` contains daily temperature, precipitation, radiation, windspeed, and CO₂.
simulation = initialize_simulation(
    cft1, initial_data;
    indices = [1],
    device = identity,       # use CuArray for a CUDA backend
    T = Float32,             # Float64 is also supported
    days = 365,
    auto_fertilizer = false,
)

run_simulation!(simulation, climate)

summary = simulation_summary(simulation)
npp = simulation.output.crop.npp
```

For a GPU simulation, construct inputs on the GPU or set `device = CuArray`
when calling `initialize_simulation`. The same process code is designed to run
over a batch of independent grid cells; `indices = [1]` selects one input grid
cell, while a longer index vector selects a larger batch.


## Design principles

### Scientific continuity, not blind replication

LPJmL provides the principal scientific baseline for the current crop-process
audit. We compare process logic, conservation behaviour, and numerical guards
where this is informative. However, exact reproduction of every LPJmL internal
variable is not the objective. Agrocosm favours documented, testable process
behaviour and makes intentional simplifications explicit.


## Testing

Run the CPU test suite with:

```bash
julia --project=. test/runtests.jl
```

GPU tests are separate because they require a functional CUDA device. For
example:

```bash
julia --project=. test/simulations/test_daily_crop_C3_precision_gpu.jl
julia --project=. test/processes/soil/test_soil_process_kernels_gpu.jl
```

## Long-term goals

When the roadmap is complete, Agrocosm should support:

- Large-domain, high-resolution crop simulation on both CPUs and GPUs
- Gradient-based calibration of cultivar and physiological parameters
- Assimilation of remotely sensed LAI, GPP, evapotranspiration, and biomass
- Combining Agrocosm with data-driven models to surpport hybrid modelling
- Sensitivity analyses for climate change and crop management
- An open software ecosystem that agricultural modellers can work together.

## Contributing

Contributions, ideas, issue reports, and process-comparison experiments are
very welcome. Agrocosm is most useful when crop physiologists, Earth system
modellers, numerical scientists, and machine-learning researchers can inspect,
question, and improve its assumptions together. Please open a GitHub issue to
start a discussion.

## Acknowledgements

Agrocosm.jl is a research project developed with the support of the
[Earth System Modeling group](https://www.asg.ed.tum.de/) at the Technical
University of Munich (TUM) and the
[FutureLab on Artificial Intelligence](https://www.pik-potsdam.de/en/futurelab-ai)
at the Potsdam Institute for Climate Impact Research (PIK). The author
acknowledges funding from the China Scholarship Council (grant agreement
no. 202303250017) and the Horizon Europe ClimTip project (grant agreement
no. 101137601).

## License

Agrocosm.jl is released under the [MIT License](LICENSE).
