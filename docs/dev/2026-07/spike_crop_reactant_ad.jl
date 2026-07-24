# Phase 6 spike: reverse-mode automatic differentiation of the crop soil biogeochemistry through
# Reactant + Enzyme. Terrarium targets full differentiability; this compiles the whole gradient of a
# scalar loss (mean-square soil nitrate after `nsteps`) with respect to the initial soil mineral-N
# state and checks it is finite and nonzero — i.e. the ported crop C–N dynamics (mineralization,
# nitrification, denitrification) differentiate end-to-end through the traced stepping loop.
#
# The differentiation target is a soil-only model carrying the crop `CropSoilBiogeochemistry`, so the
# whole rollout compiles under Reactant without the surface turbulent-flux path.
#
# Requires Julia 1.12 and a Reactant-enabled environment (Reactant, CUDA, Enzyme). Run from such an
# env, e.g.:
#   julia +1.12 --project=<ad-env> docs/dev/2026-07/spike_crop_reactant_ad.jl

using Agrocosm
using Terrarium
using Reactant, CUDA   # CUDA provides Reactant's kernel integration, even on CPU
using Enzyme
using Statistics: mean

NF = Float32

# Build the grid on `ReactantState()` so the state lives on the device and the stepping loop is traced.
grid = ColumnGrid(ReactantState(), NF, UniformSpacing(Δz = 0.1f0, N = 10))
soil = SoilEnergyWaterCarbon(NF; biogeochem = CropSoilBiogeochemistry(NF))
model = SoilModel(grid; soil)
integrator = initialize(model)

# Scalar loss: mean-square soil nitrate after `nsteps`, advanced through the traced `run_timesteps!`.
function loss(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    return mean(interior(integrator.state.soil_nitrate) .^ 2)
end

# Reverse-mode gradient; the sensitivity to the initial mineral-N pools accumulates in `dintegrator`.
function grad_loss!(integrator, dintegrator, Δt, nsteps, checkpointing)
    _, loss_value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        loss, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt), Enzyme.Const(nsteps), Enzyme.Const(checkpointing),
    )
    return loss_value
end

Δt = NF(600)
nsteps = 20
checkpointing = Reactant.Periodic(isqrt(nsteps))
dintegrator = Enzyme.make_zero(integrator)

compiled_grad! = Reactant.@compile raise = true raise_first = true sync = true grad_loss!(
    integrator, dintegrator, Δt, nsteps, checkpointing,
)
loss_value = compiled_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)

d_ammonium = Array(interior(dintegrator.state.soil_ammonium))
d_nitrate = Array(interior(dintegrator.state.soil_nitrate))

println("SPIKE OK")
println("  loss (mean nitrate²)              = ", Reactant.to_number(loss_value))
println("  max|∂loss/∂ammonium₀|             = ", maximum(abs, d_ammonium))
println("  max|∂loss/∂nitrate₀|              = ", maximum(abs, d_nitrate))

@assert all(isfinite, d_ammonium) "gradient wrt initial ammonium is finite"
@assert all(isfinite, d_nitrate) "gradient wrt initial nitrate is finite"
@assert maximum(abs, d_ammonium) > 0 || maximum(abs, d_nitrate) > 0 "the loss depends on the initial mineral-N state"
println("SPIKE ASSERTIONS PASSED")
