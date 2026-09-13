using Agrocosm
using Test

@testset "Water stress accelerates development only after flowering" begin
    accelerate = Agrocosm.drought_development_acceleration
    flowering = Float32(cft1.flowering_start)

    # CONTRACT. A zero rate is the inherited thermal-time development, whatever
    # the crop's water status or stage.
    for wscal in (0.0f0, 0.3f0, 1.0f0), fphu in (0.1f0, 0.5f0, 0.99f0)
        @test accelerate(wscal, fphu, flowering, 0.0f0) === 1.0f0
    end

    # An unstressed crop develops at the thermal-time rate at every stage.
    for fphu in (0.1f0, 0.5f0, 0.99f0)
        @test accelerate(1.0f0, fphu, flowering, 0.8f0) === 1.0f0
    end

    # Before flowering the multiplier is exactly one even under total stress.
    # Accelerating `fphu` there would run the phenological LAI curve up faster,
    # so a drought would BUILD canopy rather than shorten the season.
    @test accelerate(0.0f0, flowering - 0.01f0, flowering, 0.8f0) === 1.0f0
    @test accelerate(0.0f0, flowering, flowering, 0.8f0) > 1.0f0

    # After flowering it is APSIM's linear form, bounded by the rate.
    @test accelerate(0.5f0, 0.8f0, flowering, 0.8f0) ≈ 1.4f0
    @test accelerate(0.0f0, 0.8f0, flowering, 0.8f0) ≈ 1.8f0

    # Monotone in stress, and a sufficiency outside [0, 1] cannot push the
    # multiplier below one or past its rate.
    previous = 1.0f0
    for wscal in (1.0f0, 0.75f0, 0.5f0, 0.25f0, 0.0f0)
        current = accelerate(wscal, 0.8f0, flowering, 0.5f0)
        @test current >= previous
        previous = current
    end
    @test accelerate(1.5f0, 0.8f0, flowering, 0.5f0) === 1.0f0
    @test accelerate(-0.5f0, 0.8f0, flowering, 0.5f0) === accelerate(0.0f0, 0.8f0, flowering, 0.5f0)
end

@testset "A stressed stand reaches senescence and harvest sooner" begin
    # End to end through the kernel: same weather, same heat units, one crop
    # water-stressed after flowering and one not. The stressed stand must
    # accumulate `fphu` faster, because senescence and harvest are thresholds on
    # it and a shortened season is the whole point.
    accelerated = Agrocosm.CFTParameters{Float32, Int32}(;
        (f => (f === :drought_phenology_rate ? 1.0f0 : getfield(cft1, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)

    function advance(cft, wscal)
        crop = init_crop(1, identity)
        crop.state.phenology.is_growing .= 1
        crop.auxiliary.phenology.phu .= 2000.0f0
        crop.state.phenology.husum .= 1200.0f0   # past flowering_start
        crop.state.water.sufficiency .= wscal
        state = test_model_state(crop, init_soil(1, soilparams.soildepth, identity))
        phenology_crop!(state, Float32[60], cft, Float32[25], Float32[12])
        return crop.state.phenology.husum[1]
    end

    unstressed = advance(accelerated, 1.0f0)
    stressed = advance(accelerated, 0.2f0)
    shipped = advance(cft1, 0.2f0)

    @test stressed > unstressed
    @test unstressed === shipped
end
