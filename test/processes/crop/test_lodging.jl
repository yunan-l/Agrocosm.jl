using Agrocosm
using Test

@testset "Lodging pressure needs a storm, a filled ear and wet ground" begin
    pressure = Agrocosm.lodging_pressure_today
    at(; wind, wetness = 1.0f0, load = 300.0f0, fphu = 0.8f0,
       threshold = 8.0f0, flowering = 0.45f0, reference = 300.0f0) =
        pressure(wind, wetness, load, fphu, threshold, flowering, reference)

    # Each of the three factors can veto the day on its own.
    @test at(wind = 5.0f0) == 0.0f0                    # below the threshold
    @test at(wind = 8.0f0) == 0.0f0                    # exactly at it
    @test at(wind = 20.0f0, fphu = 0.2f0) == 0.0f0     # before flowering
    @test at(wind = 20.0f0, wetness = 0.0f0) == 0.0f0  # anchorage holds in dry ground
    @test at(wind = 20.0f0, load = 0.0f0) == 0.0f0     # nothing to overturn
    @test at(wind = 20.0f0) > 0.0f0

    # Drag goes as the square of the wind ABOVE the threshold.
    @test at(wind = 10.0f0) ≈ 4.0f0
    @test at(wind = 12.0f0) ≈ 16.0f0
    @test at(wind = 20.0f0) ≈ 144.0f0

    # Linear in load up to the reference, then saturating: a crop cannot be more
    # than fully susceptible.
    @test at(wind = 20.0f0, load = 150.0f0) ≈ 72.0f0
    @test at(wind = 20.0f0, load = 600.0f0) ≈ at(wind = 20.0f0, load = 300.0f0)
    # And linear in wetness, clamped the same way.
    @test at(wind = 20.0f0, wetness = 0.5f0) ≈ 72.0f0
    @test at(wind = 20.0f0, wetness = 1.6f0) ≈ at(wind = 20.0f0, wetness = 1.0f0)
    @test at(wind = 20.0f0, wetness = -0.3f0) == 0.0f0
    # A zero reference load cannot divide by zero; it means "always susceptible".
    @test at(wind = 20.0f0, load = 1.0f0, reference = 0.0f0) ≈ 144.0f0
end

@testset "Lodging recovery is inert at zero rate and bounded at one" begin
    recovery = Agrocosm.lodging_recovery

    # CONTRACT. Zero rate returns exactly one at any exposure, so the shipped
    # model harvests bitwise as before.
    for exposure in (0.0f0, 100.0f0, 1.0f5)
        @test recovery(exposure, 0.0f0, 0.0f0) === 1.0f0
    end

    # Below the tolerance nothing is lost, however large the rate.
    @test recovery(50.0f0, 100.0f0, 1.0f0) === 1.0f0
    # Above it, linear in the excess and clamped at total loss.
    @test recovery(600.0f0, 100.0f0, 1.0f-3) ≈ 0.5f0
    @test recovery(1.0f6, 100.0f0, 1.0f-3) === 0.0f0
    # Monotone.
    previous = 1.0f0
    for exposure in (0.0f0, 200.0f0, 400.0f0, 800.0f0, 1600.0f0)
        current = recovery(exposure, 100.0f0, 1.0f-3)
        @test current <= previous
        previous = current
    end
end

@testset "Lodging exposure accumulates over the season and resets at sowing" begin
    lodging_cft = Agrocosm.CFTParameters{Float32, Int32}(;
        (f => (f === :lodging_rate ? 1.0f-3 : getfield(cft1, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)

    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    crop.state.phenology.is_growing .= 1
    crop.auxiliary.phenology.fphu .= 0.8f0
    crop.state.carbon.storage .= 200.0f0
    crop.state.carbon.leaf .= 100.0f0
    soil.water.relative_content .= 1.0f0

    lodging!(lodging_cft, state, Float32[20.0], state)
    after_one = crop.state.phenology.lodging_exposure[1]
    @test after_one ≈ 144.0f0
    lodging!(lodging_cft, state, Float32[20.0], state)
    @test crop.state.phenology.lodging_exposure[1] ≈ 2 * after_one

    # A day outside the season contributes nothing.
    crop.state.phenology.is_growing .= 0
    lodging!(lodging_cft, state, Float32[30.0], state)
    @test crop.state.phenology.lodging_exposure[1] ≈ 2 * after_one
end
