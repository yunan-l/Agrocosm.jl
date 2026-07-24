# Reactant + Enzyme reverse-mode autodiff test for the crop soil biogeochemistry (mirrors Terrarium's
# test/reactant/autodiff.jl).
#
# Compiles a scalar loss (mean-square soil nitrate after `nsteps`) and its reverse-mode gradient with
# Reactant, and checks the sensitivity w.r.t. the initial soil mineral-nitrogen state is finite and
# nonzero — i.e. the ported crop C–N transforms (mineralization → NH₄ → nitrification → NO₃ →
# denitrification) differentiate end-to-end through the traced stepping loop. Also exercises the
# extension's `checkpointing` argument: the gradient must be identical whether the traced loop stores
# every state (`checkpointing=false`) or only periodic checkpoints — checkpointing changes memory use,
# not the math.

using Enzyme
using Statistics: mean

function ad_loss(integrator, Δt, nsteps, checkpointing)
    run_timesteps!(integrator, Δt, nsteps, checkpointing)
    return mean(interior(integrator.state.soil_nitrate) .^ 2)
end

function ad_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)
    _, loss_value = Enzyme.autodiff(
        Enzyme.set_strong_zero(Enzyme.ReverseWithPrimal),
        ad_loss, Enzyme.Active,
        Enzyme.Duplicated(integrator, dintegrator),
        Enzyme.Const(Δt),
        Enzyme.Const(nsteps),
        Enzyme.Const(checkpointing),
    )
    return loss_value
end

# Compile the gradient for a given checkpointing scheme and return
# (loss, ∂loss/∂nitrate₀ as a CPU array).
function reactant_gradient(config, NF, Δt, nsteps, checkpointing)
    integrator = build_integrator(Val(config), ReactantState(), NF)
    dintegrator = Enzyme.make_zero(integrator)
    compiled_grad! = Reactant.@compile raise = true raise_first = true sync = true ad_grad!(
        integrator, dintegrator, Δt, nsteps, checkpointing
    )
    loss_value = compiled_grad!(integrator, dintegrator, Δt, nsteps, checkpointing)
    return _scalar(loss_value), Array(interior(dintegrator.state.soil_nitrate))
end

@testset "Reactant + Enzyme autodiff (crop soil biogeochemistry)" begin
    NF = DEFAULT_NF
    Δt = NF(600)
    nsteps = 20

    @testset "gradient is finite and nonzero" begin
        loss_value, d_nitrate = reactant_gradient(:crop_soil_biogeochemistry, NF, Δt, nsteps, false)
        println("autodiff: loss=$loss_value  max|∂loss/∂NO₃₀|=$(maximum(abs, d_nitrate))")
        @test isfinite(loss_value)
        @test loss_value > 0               # mean-square nitrate is positive
        @test all(isfinite, d_nitrate)     # gradient finite everywhere
        @test maximum(abs, d_nitrate) > 0  # the loss depends on the initial nitrate state
    end

    @testset "checkpointing scheme leaves the gradient unchanged" begin
        _, d_ref = reactant_gradient(:crop_soil_biogeochemistry, NF, Δt, nsteps, false)
        _, d_ckpt = reactant_gradient(:crop_soil_biogeochemistry, NF, Δt, nsteps, Reactant.Periodic(5))
        @test all(isfinite, d_ckpt)
        @test d_ref ≈ d_ckpt               # checkpointing changes memory use, not the gradient
    end
end
