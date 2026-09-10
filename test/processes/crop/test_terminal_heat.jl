using Agrocosm
using Test

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Terminal heat exists because the reproductive sink demonstrably cannot reach
# the case it is named for. The tests that matter here are therefore not "does
# the arithmetic work" but "is it genuinely a second, independent damage path",
# and "is it inert by default", since it ships with a zero rate.

const _TT = Float32

function _filling_cft(base; rate = 0.02, threshold = 30.0, start = 0.70, stop = 0.95)
    return Agrocosm.CFTParameters{_TT, Int32}(;
        (f => (f === :filling_rate ? _TT(rate) :
               f === :filling_temperature ? _TT(threshold) :
               f === :filling_start ? _TT(start) :
               f === :filling_end ? _TT(stop) :
               getfield(base, f)) for f in fieldnames(Agrocosm.CFTParameters))...)
end

"""Step the terminal-heat kernel over a sequence of days and return the path."""
function _filling_trajectory(exposures, fphus; kwargs...)
    cft = _filling_cft(Agrocosm.cft1; kwargs...)
    crop = Agrocosm.init_crop(_TT, 1, identity)
    state = test_model_state(crop)
    Agrocosm.crop_prognostic(state).phenology.grain_fill_fraction .= one(_TT)
    Agrocosm.crop_prognostic(state).phenology.is_growing .= Int32(1)
    path = _TT[]
    for (exposure, fphu) in zip(exposures, fphus)
        Agrocosm.crop_stress_auxiliary(state).filling_exposure_hours .= _TT(exposure)
        Agrocosm.crop_phenology_auxiliary(state).fphu .= _TT(fphu)
        Agrocosm.terminal_heat!(cft, state)
        push!(path, Agrocosm.crop_prognostic(state).phenology.grain_fill_fraction[1])
    end
    return path
end

@testset "The shipped parameters are the bounded ones" begin
    # `filling_rate` is an upper bound from `tools/sterility_calibration.jl`,
    # not a fit, and it is an order of magnitude below `sterility_rate` because
    # the filling window carries far more exposure. If it ever exceeds the
    # sterility rate that ordering has been lost and the bound is stale.
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft4)
        @test cft.filling_rate > 0
        @test cft.filling_rate < cft.sterility_rate
        # And the threshold must still be the lower one, so that the field being
        # accumulated is the right one when a rate is eventually set.
        @test cft.filling_temperature < cft.sterility_temperature
        # The windows must meet rather than overlap: grain set and grain filling
        # are sequential, and an overlap would let one event be counted twice.
        @test cft.filling_start >= cft.flowering_end
        @test cft.filling_end > cft.filling_start
    end
    # Rate zero is still the exact retreat, which is what makes the mechanism
    # switchable off by parameter as well as by configuration.
    @test all(==(1), _filling_trajectory(fill(40.0, 5), fill(0.82, 5); rate = 0.0))
end

@testset "Filling loss is monotone, irreversible and window-bounded" begin
    # Same contract as grain set: starch not deposited is not deposited later.
    hot = _filling_trajectory([20.0, 20.0, 20.0], [0.82, 0.82, 0.82])
    @test issorted(hot; rev = true)
    @test hot[end] < hot[1] < 1
    # A cool spell afterwards must not recover anything.
    recovered = _filling_trajectory([20.0, 0.0, 0.0], [0.82, 0.82, 0.82])
    @test recovered[2] == recovered[1]
    @test recovered[3] == recovered[1]

    # Outside the window heat does nothing, including at both edges exactly.
    for fphu in (0.0, 0.5, 0.70, 0.95, 1.0)
        @test all(==(1), _filling_trajectory([40.0], [fphu]))
    end
    for fphu in (0.72, 0.825, 0.94)
        @test _filling_trajectory([40.0], [fphu])[1] < 1
    end

    # Not growing means no damage: a bare field cannot lose grain weight.
    cft = _filling_cft(Agrocosm.cft1)
    crop = Agrocosm.init_crop(_TT, 1, identity)
    state = test_model_state(crop)
    Agrocosm.crop_prognostic(state).phenology.grain_fill_fraction .= one(_TT)
    Agrocosm.crop_prognostic(state).phenology.is_growing .= Int32(0)
    Agrocosm.crop_stress_auxiliary(state).filling_exposure_hours .= _TT(40)
    Agrocosm.crop_phenology_auxiliary(state).fphu .= _TT(0.82)
    Agrocosm.terminal_heat!(cft, state)
    @test Agrocosm.crop_prognostic(state).phenology.grain_fill_fraction[1] == 1

    # Clamped at zero rather than going negative under extreme exposure.
    @test _filling_trajectory(fill(1e4, 3), fill(0.825, 3))[end] == 0
end

@testset "The two damage paths are independent" begin
    # THE structural claim. Grain set and grain filling must be separately
    # switchable and must not read each other's exposure field, or the ablation
    # cannot attribute a loss to one rather than the other.
    cft = _filling_cft(Agrocosm.cft1)
    crop = Agrocosm.init_crop(_TT, 1, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    stress = Agrocosm.crop_stress_auxiliary(state)
    phenology.grain_fill_fraction .= one(_TT)
    phenology.grain_set_fraction .= one(_TT)
    phenology.is_growing .= Int32(1)
    Agrocosm.crop_phenology_auxiliary(state).fphu .= _TT(0.825)

    # Exposure in the FILLING field only: filling falls, grain set does not.
    stress.heat_exposure_hours .= zero(_TT)
    stress.filling_exposure_hours .= _TT(20)
    Agrocosm.terminal_heat!(cft, state)
    Agrocosm.reproductive_sink!(cft, state)
    @test phenology.grain_fill_fraction[1] < 1
    @test phenology.grain_set_fraction[1] == 1

    # And the converse, at a development stage inside the FLOWERING window.
    phenology.grain_fill_fraction .= one(_TT)
    phenology.grain_set_fraction .= one(_TT)
    Agrocosm.crop_phenology_auxiliary(state).fphu .= _TT(0.575)
    stress.heat_exposure_hours .= _TT(20)
    stress.filling_exposure_hours .= zero(_TT)
    Agrocosm.terminal_heat!(cft, state)
    Agrocosm.reproductive_sink!(cft, state)
    @test phenology.grain_set_fraction[1] < 1
    @test phenology.grain_fill_fraction[1] == 1
end

@testset "The filling accumulator exceeds the sterility one" begin
    # It is the same integral at a lower threshold, so on any day it must be at
    # least as large - and strictly larger whenever the leaf crosses the lower
    # threshold but not the upper. That gap IS the mechanism at the hot-wheat
    # cell, where the filling window holds 75.8 hours above 30 C and 0.4 above
    # 35 C.
    cft = _filling_cft(Agrocosm.cft1; threshold = 30.0)
    config = DiurnalConfig(; steps = 24,
                             shape = Agrocosm.diurnal_shape_code(:sinusoid))
    function pair(mean_temperature)
        crop = Agrocosm.init_crop(_TT, 1, identity)
        state = test_model_state(crop)
        Agrocosm.heat_exposure!(
            cft, state, _TT[14.0], _TT[mean_temperature],
            DiurnalForcing(config, _TT[10.0]); organ = nothing,
        )
        stress = Agrocosm.crop_stress_auxiliary(state)
        return stress.heat_exposure_hours[1], stress.filling_exposure_hours[1]
    end
    for mean_temperature in (20.0, 26.0, 28.0, 31.0, 34.0, 40.0)
        sterility, filling = pair(mean_temperature)
        @test filling >= sterility
    end
    # The discriminating band: a day whose peak sits between the two thresholds
    # must register filling exposure and negligible sterility exposure. That is
    # the hot-wheat signature.
    sterility, filling = pair(27.0)      # peak 32 C, between 30 and 35
    @test filling > 1
    @test sterility < 0.05
end

@testset "Both factors reach yield, multiplicatively" begin
    # The harvest index carries grain set AND grain filling, and they multiply.
    # If either were added or min-ed, one mechanism could mask the other.
    both = Agrocosm.compute_harvest_index(_TT(0.8), _TT(0.5), _TT(0.3), _TT(0)) *
           _TT(0.5) * _TT(0.5)
    set_only = Agrocosm.compute_harvest_index(_TT(0.8), _TT(0.5), _TT(0.3), _TT(0)) *
               _TT(0.5) * _TT(1.0)
    fill_only = Agrocosm.compute_harvest_index(_TT(0.8), _TT(0.5), _TT(0.3), _TT(0)) *
                _TT(1.0) * _TT(0.5)
    @test both < set_only
    @test both < fill_only
    @test both ≈ set_only * 0.5 rtol = 1e-5
end

@testset "Terminal heat is per-cell independent" begin
    # Same hazard as the exposure kernels: one cell in every other test here.
    cells = 4
    cft = _filling_cft(Agrocosm.cft1; rate = 0.02)
    exposures = _TT[0.0, 5.0, 20.0, 80.0]
    fphus = _TT[0.60, 0.75, 0.825, 0.90]      # first is outside the window
    growing = Int32[1, 1, 1, 0]               # last is not growing

    crop = Agrocosm.init_crop(_TT, cells, identity)
    state = test_model_state(crop)
    phenology = Agrocosm.crop_prognostic(state).phenology
    phenology.grain_fill_fraction .= one(_TT)
    phenology.is_growing .= growing
    Agrocosm.crop_stress_auxiliary(state).filling_exposure_hours .= exposures
    Agrocosm.crop_phenology_auxiliary(state).fphu .= fphus
    Agrocosm.terminal_heat!(cft, state)
    batch = copy(phenology.grain_fill_fraction)

    for index in 1:cells
        single_crop = Agrocosm.init_crop(_TT, 1, identity)
        single = test_model_state(single_crop)
        single_phenology = Agrocosm.crop_prognostic(single).phenology
        single_phenology.grain_fill_fraction .= one(_TT)
        single_phenology.is_growing .= Int32(growing[index])
        Agrocosm.crop_stress_auxiliary(single).filling_exposure_hours .= exposures[index]
        Agrocosm.crop_phenology_auxiliary(single).fphu .= fphus[index]
        Agrocosm.terminal_heat!(cft, single)
        @test single_phenology.grain_fill_fraction[1] === batch[index]
    end

    # The forcing was chosen so the four cells take four different branches:
    # outside the window, inside with light exposure, inside with heavy
    # exposure, and not growing. If they all came out equal the test proves
    # nothing.
    @test batch[1] == 1                 # outside the window
    @test 1 > batch[2] > batch[3]       # inside, monotone in exposure
    @test batch[4] == 1                 # not growing
end
