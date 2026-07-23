# Phase 0 toolchain spike: run a Terrarium SoilModel column on CPU.
#
# Purpose: validate that the Terrarium framework (grid -> model -> initialize ->
# timestep!/run!) works end-to-end on this machine before Agrocosm's own
# infrastructure is deleted in Phase 1.
#
# Run from the Terrarium.jl project environment:
#   julia --project=/Users/max/Nextcloud/Terrarium.jl spike_soil_column_cpu.jl
#
# GPU note: a GPU spike (ColumnGrid(GPU(), ...)) cannot be executed on this
# machine (Apple Silicon, no CUDA). GPU validation is deferred to CI/a CUDA host.

using Terrarium

const NF = Float32

# 1D column with 10 exponentially spaced soil layers on the CPU.
grid = ColumnGrid(CPU(), NF, ExponentialSpacing(N = 10))
@assert eltype(grid) === NF

# Quasi-steady-state temperature init + fully saturated water/ice.
initializer = SoilInitializer(
    eltype(grid),
    energy = QuasiThermalSteadyState(eltype(grid), T₀ = -1.0),
    hydrology = ConstantSaturation(eltype(grid), sat = 1.0),
)

model = SoilModel(grid; timestepper = ForwardEuler(eltype(grid)), initializer = initializer)

boundary_conditions = PrescribedSurfaceTemperature(:T_ub, 1.0)

integrator = initialize(model; boundary_conditions)

# One warm-up step, then a short integration.
timestep!(integrator)
run!(integrator, period = Day(3))

T = interior(integrator.state.temperature)[1, 1, :]
f = interior(integrator.state.liquid_water_fraction)[1, 1, :]
zs = znodes(integrator.state.temperature)

println("SPIKE OK")
println("  eltype(grid)      = ", eltype(grid))
println("  n layers          = ", length(T))
println("  current_time      = ", current_time(integrator))
println("  temperature[1:3]  = ", T[1:min(3, end)])
println("  liquid_frac[1:3]  = ", f[1:min(3, end)])
println("  z nodes[1:3]      = ", zs[1:min(3, end)])

@assert all(isfinite, T) "non-finite temperatures"
@assert all(isfinite, f) "non-finite liquid fractions"
@assert length(T) == 10
println("SPIKE ASSERTIONS PASSED")
