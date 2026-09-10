using Agrocosm
using Test

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# The claims worth testing here are structural, not arithmetic: that water and
# heat share one floret pool, that the order they are applied in cannot matter,
# and that the mechanism reads the DAILY sufficiency rather than the
# season-cumulative deficit field whose sentinel value would read a bare field
# as total drought.

const _WT = Float32

function _water_cft(base; rate = 0.05, sufficiency = 0.5)
    return Agrocosm.CFTParameters{_WT, Int32}(;
        (f => (f === :water_sterility_rate ? _WT(rate) :
               f === :water_sterility_sufficiency ? _WT(sufficiency) :
               getfield(base, f)) for f in fieldnames(Agrocosm.CFTParameters))...)
end

"""State with one cell in the flowering window and a given daily sufficiency."""
function _stand(cft; sufficiency = 0.2, fphu = 0.575, growing = true,
                exposure = 0.0, set = 1.0)
    crop = Agrocosm.init_crop(_WT, 1, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    phenology.grain_set_fraction .= _WT(set)
    phenology.is_growing .= Int32(growing)
    Agrocosm.crop_prognostic(state).water.sufficiency .= _WT(sufficiency)
    Agrocosm.crop_phenology_auxiliary(state).fphu .= _WT(fphu)
    Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours .= _WT(exposure)
    return state
end

@testset "Water sterility ships inert and needs no prerequisite" begin
    # cft9 is soybean, which the five-cell gate actually runs; cft4 is a
    # crop it never runs. The defaults are shared by every CFT, so this
    # is about sampling the ones the project uses.
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.water_sterility_rate == 0
        @test 0 < cft.water_sterility_sufficiency <= 1
    end
    # No exposure source to enable: unlike the heat mechanisms, this reads a
    # field every configuration writes. A configuration with nothing else on
    # must build.
    configuration = Agrocosm.SimulationConfiguration(
        Float32, identity, 10, [1], [1]; water_sterility = true,
    )
    @test configuration.water_sterility
    @test Agrocosm.diurnal_configuration(configuration) === nothing
    @test Agrocosm.heat_exposure_configuration(configuration) === nothing

    state = _stand(_water_cft(Agrocosm.cft1; rate = 0.0))
    Agrocosm.water_sterility!(_water_cft(Agrocosm.cft1; rate = 0.0), state)
    @test Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1] == 1
end

@testset "Loss is driven by the daily shortfall, in the flowering window" begin
    cft = _water_cft(Agrocosm.cft1)
    # Above the threshold: nothing, exactly.
    for sufficiency in (0.5, 0.7, 1.0)
        state = _stand(cft; sufficiency)
        Agrocosm.water_sterility!(cft, state)
        @test Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1] == 1
    end
    # Below it: monotone in the shortfall.
    losses = map((0.4, 0.3, 0.2, 0.0)) do sufficiency
        state = _stand(cft; sufficiency)
        Agrocosm.water_sterility!(cft, state)
        return 1 - Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1]
    end
    @test issorted(losses)
    @test all(>(0), losses)

    # Outside the window, and on a bare field, nothing happens - the latter
    # because `water.sufficiency` is 1 for an absent stand, and because the
    # kernel gates on is_growing regardless.
    for fphu in (0.0, 0.45, 0.70, 1.0)
        state = _stand(cft; sufficiency = 0.0, fphu)
        Agrocosm.water_sterility!(cft, state)
        @test Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1] == 1
    end
    state = _stand(cft; sufficiency = 0.0, growing = false)
    Agrocosm.water_sterility!(cft, state)
    @test Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1] == 1

    # Monotone and irreversible across days.
    cft2 = _water_cft(Agrocosm.cft1)
    state = _stand(cft2; sufficiency = 0.2)
    path = _WT[]
    for sufficiency in (0.2, 0.2, 1.0, 1.0)
        Agrocosm.crop_prognostic(state).water.sufficiency .= _WT(sufficiency)
        Agrocosm.water_sterility!(cft2, state)
        push!(path, Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1])
    end
    @test issorted(path; rev = true)
    @test path[3] == path[2]        # the drought broke; nothing recovers
    @test path[4] == path[2]
    @test path[2] < path[1] < 1
end

@testset "Water and heat share one floret pool, order-independently" begin
    # THE structural claim. Both subtract from `grain_set_fraction`, so the
    # order cannot matter and the two losses cannot exceed the pool.
    cft = _water_cft(Agrocosm.cft1)
    hot = Agrocosm.CFTParameters{_WT, Int32}(;
        (f => (f === :sterility_rate ? _WT(0.02) :
               f === :water_sterility_rate ? _WT(0.05) :
               f === :water_sterility_sufficiency ? _WT(0.5) :
               getfield(Agrocosm.cft1, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)

    water_first = _stand(hot; sufficiency = 0.2, exposure = 6.0)
    Agrocosm.water_sterility!(hot, water_first)
    Agrocosm.reproductive_sink!(hot, water_first)
    heat_first = _stand(hot; sufficiency = 0.2, exposure = 6.0)
    Agrocosm.reproductive_sink!(hot, heat_first)
    Agrocosm.water_sterility!(hot, heat_first)
    @test Agrocosm.crop_prognostic(water_first).phenology.grain_set_fraction[1] ===
          Agrocosm.crop_prognostic(heat_first).phenology.grain_set_fraction[1]

    # Both together must cost more than either alone, and never go below zero.
    water_only = _stand(hot; sufficiency = 0.2, exposure = 0.0)
    Agrocosm.water_sterility!(hot, water_only)
    Agrocosm.reproductive_sink!(hot, water_only)
    heat_only = _stand(hot; sufficiency = 1.0, exposure = 6.0)
    Agrocosm.water_sterility!(hot, heat_only)
    Agrocosm.reproductive_sink!(hot, heat_only)
    both = Agrocosm.crop_prognostic(water_first).phenology.grain_set_fraction[1]
    @test both < Agrocosm.crop_prognostic(water_only).phenology.grain_set_fraction[1]
    @test both < Agrocosm.crop_prognostic(heat_only).phenology.grain_set_fraction[1]

    # Extreme of both: clamped at zero, not negative, in either order.
    ruin = Agrocosm.CFTParameters{_WT, Int32}(;
        (f => (f === :sterility_rate ? _WT(10.0) :
               f === :water_sterility_rate ? _WT(10.0) :
               f === :water_sterility_sufficiency ? _WT(0.5) :
               getfield(Agrocosm.cft1, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
    for order in ((Agrocosm.water_sterility!, Agrocosm.reproductive_sink!),
                  (Agrocosm.reproductive_sink!, Agrocosm.water_sterility!))
        state = _stand(ruin; sufficiency = 0.0, exposure = 100.0)
        order[1](ruin, state)
        order[2](ruin, state)
        @test Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1] == 0
    end
end

@testset "Water sterility is per-cell independent" begin
    cells = 4
    cft = _water_cft(Agrocosm.cft1)
    sufficiencies = _WT[1.0, 0.45, 0.2, 0.0]
    fphus = _WT[0.575, 0.575, 0.575, 0.30]      # last is outside the window

    crop = Agrocosm.init_crop(_WT, cells, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    phenology.grain_set_fraction .= one(_WT)
    phenology.is_growing .= Int32(1)
    Agrocosm.crop_prognostic(state).water.sufficiency .= sufficiencies
    Agrocosm.crop_phenology_auxiliary(state).fphu .= fphus
    Agrocosm.water_sterility!(cft, state)
    batch = copy(phenology.grain_set_fraction)

    for index in 1:cells
        single = _stand(cft; sufficiency = sufficiencies[index], fphu = fphus[index])
        Agrocosm.water_sterility!(cft, single)
        @test Agrocosm.crop_prognostic(single).phenology.grain_set_fraction[1] ===
              batch[index]
    end
    @test batch[1] == 1                 # above the threshold
    @test 1 > batch[2] > batch[3]       # monotone in the shortfall
    @test batch[4] == 1                 # outside the window
end
