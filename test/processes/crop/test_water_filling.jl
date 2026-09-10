using Agrocosm
using Test

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# The fourth cell of the {heat, water} x {number, weight} square. What is worth
# testing is the square's closure properties: the four mechanisms pair up onto
# two shared states, each pair is order-independent, and each channel touches
# only its own state and window.

const _FT = Float32

function _wf_cft(base; rate = 0.05, sufficiency = 0.5, heat_fill = 0.0,
                 water_set = 0.0, heat_set = 0.0)
    return Agrocosm.CFTParameters{_FT, Int32}(;
        (f => (f === :water_filling_rate ? _FT(rate) :
               f === :water_filling_sufficiency ? _FT(sufficiency) :
               f === :filling_rate ? _FT(heat_fill) :
               f === :water_sterility_rate ? _FT(water_set) :
               f === :sterility_rate ? _FT(heat_set) :
               getfield(base, f)) for f in fieldnames(Agrocosm.CFTParameters))...)
end

function _stand(; sufficiency = 0.2, fphu = 0.825, growing = true,
                filling_exposure = 0.0, heat_exposure = 0.0)
    crop = Agrocosm.init_crop(_FT, 1, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    phenology.grain_set_fraction .= one(_FT)
    phenology.grain_fill_fraction .= one(_FT)
    phenology.is_growing .= Int32(growing)
    Agrocosm.crop_prognostic(state).water.sufficiency .= _FT(sufficiency)
    Agrocosm.crop_phenology_auxiliary(state).fphu .= _FT(fphu)
    stress = Agrocosm.crop_stress_auxiliary(state)
    stress.filling_exposure_hours .= _FT(filling_exposure)
    stress.heat_exposure_hours .= _FT(heat_exposure)
    return state
end

_fill(state) = Agrocosm.crop_prognostic(state).phenology.grain_fill_fraction[1]
_set(state) = Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1]

@testset "Water filling ships inert and needs no prerequisite" begin
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft4)
        @test cft.water_filling_rate == 0
        @test 0 < cft.water_filling_sufficiency <= 1
    end
    configuration = Agrocosm.SimulationConfiguration(
        Float32, identity, 10, [1], [1]; water_filling = true,
    )
    @test configuration.water_filling
    state = _stand()
    Agrocosm.water_filling!(_wf_cft(Agrocosm.cft1; rate = 0.0), state)
    @test _fill(state) == 1
end

@testset "Loss follows the daily shortfall inside the filling window" begin
    cft = _wf_cft(Agrocosm.cft1)
    for sufficiency in (0.5, 0.8, 1.0)
        state = _stand(; sufficiency)
        Agrocosm.water_filling!(cft, state)
        @test _fill(state) == 1
    end
    losses = map((0.4, 0.3, 0.1, 0.0)) do sufficiency
        state = _stand(; sufficiency)
        Agrocosm.water_filling!(cft, state)
        return 1 - _fill(state)
    end
    @test issorted(losses)
    @test all(>(0), losses)

    # The FILLING window, not the flowering one: a dry anthesis must not reach
    # grain weight through this channel.
    for fphu in (0.0, 0.575, 0.70, 0.95, 1.0)
        state = _stand(; sufficiency = 0.0, fphu)
        Agrocosm.water_filling!(cft, state)
        @test _fill(state) == 1
    end
    for fphu in (0.72, 0.825, 0.94)
        state = _stand(; sufficiency = 0.0, fphu)
        Agrocosm.water_filling!(cft, state)
        @test _fill(state) < 1
    end

    # Not growing, and irreversible once the drought breaks.
    state = _stand(; sufficiency = 0.0, growing = false)
    Agrocosm.water_filling!(cft, state)
    @test _fill(state) == 1

    state = _stand(; sufficiency = 0.2)
    path = _FT[]
    for sufficiency in (0.2, 0.2, 1.0, 1.0)
        Agrocosm.crop_prognostic(state).water.sufficiency .= _FT(sufficiency)
        Agrocosm.water_filling!(cft, state)
        push!(path, _fill(state))
    end
    @test issorted(path; rev = true)
    @test path[3] == path[2] == path[4]
    @test path[2] < path[1] < 1
end

@testset "The 2x2 square closes onto two states" begin
    # Each channel must touch ONLY its own state, so the ablation can attribute
    # a loss to grain number or grain weight without ambiguity.
    all_on = _wf_cft(Agrocosm.cft1; rate = 0.05, heat_fill = 0.02,
                     water_set = 0.05, heat_set = 0.02)

    # In the FILLING window with both drivers live: only grain weight moves.
    state = _stand(; fphu = 0.825, sufficiency = 0.2, filling_exposure = 10.0,
                   heat_exposure = 10.0)
    Agrocosm.water_filling!(all_on, state)
    Agrocosm.terminal_heat!(all_on, state)
    Agrocosm.water_sterility!(all_on, state)
    Agrocosm.reproductive_sink!(all_on, state)
    @test _fill(state) < 1
    @test _set(state) == 1

    # In the FLOWERING window: only grain number moves.
    state = _stand(; fphu = 0.575, sufficiency = 0.2, filling_exposure = 10.0,
                   heat_exposure = 10.0)
    Agrocosm.water_filling!(all_on, state)
    Agrocosm.terminal_heat!(all_on, state)
    Agrocosm.water_sterility!(all_on, state)
    Agrocosm.reproductive_sink!(all_on, state)
    @test _set(state) < 1
    @test _fill(state) == 1
end

@testset "Water and heat filling are order-independent" begin
    cft = _wf_cft(Agrocosm.cft1; rate = 0.05, heat_fill = 0.02)
    water_first = _stand(; sufficiency = 0.2, filling_exposure = 10.0)
    Agrocosm.water_filling!(cft, water_first)
    Agrocosm.terminal_heat!(cft, water_first)
    heat_first = _stand(; sufficiency = 0.2, filling_exposure = 10.0)
    Agrocosm.terminal_heat!(cft, heat_first)
    Agrocosm.water_filling!(cft, heat_first)
    @test _fill(water_first) === _fill(heat_first)

    # Both cost more than either alone, and the pool clamps at zero.
    water_only = _stand(; sufficiency = 0.2, filling_exposure = 0.0)
    Agrocosm.water_filling!(cft, water_only)
    Agrocosm.terminal_heat!(cft, water_only)
    heat_only = _stand(; sufficiency = 1.0, filling_exposure = 10.0)
    Agrocosm.water_filling!(cft, heat_only)
    Agrocosm.terminal_heat!(cft, heat_only)
    @test _fill(water_first) < _fill(water_only)
    @test _fill(water_first) < _fill(heat_only)

    ruin = _wf_cft(Agrocosm.cft1; rate = 10.0, heat_fill = 10.0)
    for order in ((Agrocosm.water_filling!, Agrocosm.terminal_heat!),
                  (Agrocosm.terminal_heat!, Agrocosm.water_filling!))
        state = _stand(; sufficiency = 0.0, filling_exposure = 100.0)
        order[1](ruin, state); order[2](ruin, state)
        @test _fill(state) == 0
    end
end

@testset "Water filling is per-cell independent" begin
    cells = 4
    cft = _wf_cft(Agrocosm.cft1)
    sufficiencies = _FT[1.0, 0.45, 0.1, 0.0]
    fphus = _FT[0.825, 0.825, 0.825, 0.575]     # last is outside the window

    crop = Agrocosm.init_crop(_FT, cells, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    phenology.grain_fill_fraction .= one(_FT)
    phenology.is_growing .= Int32(1)
    Agrocosm.crop_prognostic(state).water.sufficiency .= sufficiencies
    Agrocosm.crop_phenology_auxiliary(state).fphu .= fphus
    Agrocosm.water_filling!(cft, state)
    batch = copy(phenology.grain_fill_fraction)

    for index in 1:cells
        single = _stand(; sufficiency = sufficiencies[index], fphu = fphus[index])
        Agrocosm.water_filling!(cft, single)
        @test _fill(single) === batch[index]
    end
    @test batch[1] == 1
    @test 1 > batch[2] > batch[3]
    @test batch[4] == 1
end
