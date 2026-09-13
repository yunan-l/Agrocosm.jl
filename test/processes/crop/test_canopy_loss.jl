using Agrocosm
using Test

@testset "A stressed canopy sheds leaf area, and a zero rate sheds none" begin
    shed = Agrocosm.stress_canopy_loss

    # CONTRACT. A zero rate is the inherited canopy at every stress level.
    for stress in (0.0f0, 0.3f0, 0.7f0, 1.0f0)
        @test shed(stress, 0.0f0) === 0.0f0
    end

    # An unstressed crop sheds nothing at any rate.
    @test shed(1.0f0, 0.5f0) === 0.0f0

    # Linear in the deficit, and monotone.
    @test shed(0.5f0, 0.08f0) ≈ 0.04f0
    @test shed(0.0f0, 0.08f0) ≈ 0.08f0
    previous = 0.0f0
    for stress in (1.0f0, 0.75f0, 0.5f0, 0.25f0, 0.0f0)
        current = shed(stress, 0.1f0)
        @test current >= previous
        previous = current
    end

    # A sufficiency outside [0, 1] cannot shed a negative fraction or more than
    # the whole canopy, and neither can an unbounded rate.
    @test shed(1.4f0, 0.1f0) === 0.0f0
    @test shed(-0.2f0, 0.1f0) === shed(0.0f0, 0.1f0)
    @test shed(0.0f0, 3.0f0) === 1.0f0
end

@testset "Shed leaf area is not regrown tomorrow" begin
    # The loss leaves the standing canopy while `lai_previous_potential` keeps
    # tracking the phenological trajectory, so a stressed stand stays smaller
    # than an unstressed one rather than catching up.
    with_rate(rate) = Agrocosm.CFTParameters{Float32, Int32}(;
        (f => (f === :stress_canopy_loss_rate ? Float32(rate) : getfield(cft1, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)

    function canopy_after(cft, wscal, days)
        crop = init_crop(1, identity)
        crop.state.phenology.is_growing .= 1
        crop.state.water.sufficiency .= Float32(wscal)
        crop.state.nitrogen.sufficiency .= 1.0f0
        crop.auxiliary.canopy.flaimax .= 0.5f0
        crop.state.canopy.lai .= 2.0f0
        crop.state.canopy.lai_previous_potential .= 2.0f0
        state = test_model_state(crop, init_soil(1, soilparams.soildepth, identity))
        for _ in 1:days
            lai_crop!(state, cft)
        end
        return crop.state.canopy.lai[1]
    end

    # Stress still slows EXPANSION at rate zero - that is inherited - so the
    # contract is against the shipped parameters, not against an unstressed run.
    shipped = canopy_after(cft1, 0.3f0, 10)
    @test canopy_after(with_rate(0.0), 0.3f0, 10) === shipped
    stressed = canopy_after(with_rate(0.05), 0.3f0, 10)
    @test stressed < shipped
    # Ten more days of shedding compounds rather than resetting.
    @test canopy_after(with_rate(0.05), 0.3f0, 20) < stressed
end
