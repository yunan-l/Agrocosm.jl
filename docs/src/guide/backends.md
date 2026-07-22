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

CUDA is loaded by the package but a working NVIDIA device is required only
when constructing GPU arrays or running GPU tests. Always validate a new GPU,
driver, Julia, or CUDA.jl combination with the end-to-end CPU/GPU equivalence
test before production runs.
