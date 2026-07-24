# # Differentiating the crop soil biogeochemistry with Enzyme + Reactant
#
# This combines the previous two examples: we differentiate a full model rollout (as in
# `differentiating_crop_soil.jl`) but through the Reactant-compiled program (as in
# `crop_soil_reactant.jl`). Enzyme differentiates the compiled StableHLO rather than Julia IR, and the
# whole gradient — trace, forward pass, and reverse pass — is compiled in one go with
# `@compile`. Reactant's own `@trace checkpointing=true` loop bounds the reverse-pass memory (selected
# with the `checkpointing` argument to `run_timesteps!`), so we do not need Checkpointing.jl here.
#
# Run from this environment:  julia --project=examples/autodiff examples/autodiff/differentiating_crop_soil_reactant.jl

using Agrocosm, Terrarium
using Reactant, CUDA   # CUDA provides Reactant's kernel integration, even on CPU
using Enzyme
using Statistics: mean

NF = Float32

# The same soil column carrying the crop C–N biogeochemistry, on `ReactantState()`.
grid = ColumnGrid(ReactantState(), NF, UniformSpacing(Δz = NF(0.1), N = 10))
soil = SoilEnergyWaterCarbon(NF; biogeochem = CropSoilBiogeochemistry(NF))
model = SoilModel(grid; soil)
bcs = PrescribedSurfaceTemperature(:T_ub, NF(15))
initializers = (temperature = (x, z) -> NF(15) - NF(0.02) * z,)
integrator = initialize(model; boundary_conditions = bcs, initializers)

# Scalar objective: the mean-square soil nitrate after `nsteps`, advanced through the traced
# `run_timesteps!` (the stepping loop underlying `run!`; on `ReactantState` it takes a checkpointing
# scheme).
function loss(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    return mean(interior(integrator.state.soil_nitrate) .^ 2)
end

# Reverse mode computes the objective together with its gradient; passing the integrator as
# `Duplicated` accumulates the sensitivity of the loss w.r.t. the initial state into `dintegrator`.
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
checkpointing = Reactant.Periodic(isqrt(nsteps))   # ≈ √n checkpoints
dintegrator = Enzyme.make_zero(integrator)

# Compile the whole gradient in one program, then run it.
compiled_grad! = @compile raise = true raise_first = true sync = true grad_loss!(
    integrator, dintegrator, Δt, nsteps, checkpointing,
)
loss_value = compiled_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)

# The sensitivity of the loss w.r.t. the initial soil nitrate, by depth.
d_nitrate = Array(interior(dintegrator.state.soil_nitrate))[1, 1, :]

println("loss (mean nitrate²)                = ", Reactant.to_number(loss_value))
println("∂loss/∂(initial nitrate) by depth   = ", round.(d_nitrate, sigdigits = 3))
