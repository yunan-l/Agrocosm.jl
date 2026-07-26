# CPU, GPU, and precision

Backend arrays are selected during initialization:

```julia
# CPU
simulation = initialize_simulation(cft1, data;
    indices = [1], device = identity, T = Float64, days = 365)

# NVIDIA GPU
using CUDA
simulation_gpu = initialize_simulation(cft1, data;
    indices = 1:1000, device = CuArray, T = Float32, days = 365)
```

The same daily process code launches backend-neutral kernels through
KernelAbstractions.jl. `Float32` is generally preferable for GPU throughput;
`Float64` is useful for numerical audits and supported hardware.

Grid cells are independent batch members. Increasing the selected `indices`
expands the batch; it does not add lateral exchange among cells.

CUDA is an optional dependency: CPU users do not load it with Agrocosm.
Install and load CUDA.jl only when constructing CUDA arrays. Always validate a
new accelerator backend, driver, Julia, or package combination with the
end-to-end backend-equivalence test before production runs.

Every initialized simulation also exposes a typed execution contract:

```julia
simulation.config.execution.domain
architecture_name(simulation.config.execution)
float_type(simulation.config.execution)
```

The domain contains the stable active-cell order used by state, forcing,
output, and checkpoints.
