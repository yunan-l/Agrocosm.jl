using Agrocosm
using Test

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# `heat_exposure!` exists to break a coupling, not to add a process: it fills
# `heat_exposure_hours` without running a sub-daily assimilation loop, so the
# reproductive sink can ride on the calibrated daily kernel. The claim that
# makes it legitimate is that it is the SAME integral the sub-daily kernels
# accumulate inline - if the two ever disagree, one of them is wrong and every
# comparison between the two architectures is void. That equality is the first
# testset here and it is the reason this file exists; the rest are the ordinary
# boundary cases.

const _T = Float32

"""CFT with a sterility threshold low enough that the window actually bites."""
function _exposure_cft(base; threshold = 20.0)
    return Agrocosm.CFTParameters{_T, Int32}(;
        (f => (f === :sterility_temperature ? _T(threshold) : getfield(base, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
end

_organ(lai) = Agrocosm.OrganTemperatureForcing(
    _T[0.006], _T[101325.0], _T[2.0], _T[250.0], _T[-60.0], _T[0.2], _T[lai], _T[3.0],
)

"""Exposure written by the standalone pass for one day of forcing."""
function _standalone(cft, config, range, daylength, temperature; organ = nothing)
    crop = Agrocosm.init_crop(_T, 1, identity)
    state = test_model_state(crop)
    stress = Agrocosm.crop_stress_auxiliary(state)
    stress.heat_exposure_hours .= _T(-1)
    stress.filling_exposure_hours .= _T(-1)
    Agrocosm.heat_exposure!(
        cft, state, _T[daylength], _T[temperature],
        DiurnalForcing(config, _T[range]); organ,
    )
    return stress.heat_exposure_hours[1]
end

"""Both accumulators from the standalone pass, as a pair."""
function _standalone_pair(cft, config, range, daylength, temperature; organ = nothing)
    crop = Agrocosm.init_crop(_T, 1, identity)
    state = test_model_state(crop)
    stress = Agrocosm.crop_stress_auxiliary(state)
    stress.heat_exposure_hours .= _T(-1)
    stress.filling_exposure_hours .= _T(-1)
    Agrocosm.heat_exposure!(
        cft, state, _T[daylength], _T[temperature],
        DiurnalForcing(config, _T[range]); organ,
    )
    return stress.heat_exposure_hours[1], stress.filling_exposure_hours[1]
end

"""Exposure written by the sub-daily assimilation kernel for the same day."""
function _from_kernel(pathway, cft, config, range, daylength, temperature;
                      organ = nothing, apar = 5.0, stress = 0.8)
    return first(_kernel_pair(pathway, cft, config, range, daylength, temperature;
                              organ, apar, stress))
end

"""Both accumulators from the sub-daily assimilation kernel, as a pair."""
function _kernel_pair(pathway, cft, config, range, daylength, temperature;
                      organ = nothing, apar = 5.0, stress = 0.8)
    crop = Agrocosm.init_crop(_T, 1, identity)
    state = test_model_state(crop)
    crop.auxiliary.photosynthesis.temperature_stress .= _T(stress)
    aux = Agrocosm.crop_stress_auxiliary(state)
    aux.heat_exposure_hours .= _T(-1)
    aux.filling_exposure_hours .= _T(-1)
    photosynthesis!(
        pathway, cft, state, _T[apar], _T[daylength], _T[temperature], _T[400.0],
        DiurnalForcing(config, _T[range]); comp_vcmax = true, organ,
    )
    return aux.heat_exposure_hours[1], aux.filling_exposure_hours[1]
end

@testset "The standalone pass reproduces the sub-daily kernel's integral" begin
    # THE drift guard. Both routes are given the same daily state, so they must
    # produce the same number - bitwise, because it is the same expression
    # accumulated in the same order over the same trip count, not an
    # approximation of it. Anything less than exact equality here means the two
    # architectures are no longer comparable and the exposure-based ablation
    # cells stop measuring what they claim.
    for steps in (1, 12, 24, 48)
        for shape in (:flat, :sinusoid, :daytime_neutral)
            config = DiurnalConfig(; steps,
                                     shape = Agrocosm.diurnal_shape_code(shape))
            for (pathway, base) in ((Val(:C3), Agrocosm.cft1), (Val(:C4), Agrocosm.cft3))
                cft = _exposure_cft(base)
                for (range, daylength, temperature) in (
                    (12.0, 14.0, 24.0),   # threshold crossed for part of the day
                    (0.0, 12.0, 30.0),    # no range: every sub-step identical
                    (20.0, 10.0, 15.0),   # cold mean, wide range
                    (6.0, 16.0, 40.0),    # hot all day
                )
                    for organ in (nothing, _organ(3.0))
                        standalone = _standalone_pair(cft, config, range, daylength,
                                                      temperature; organ)
                        kernel = _kernel_pair(pathway, cft, config, range, daylength,
                                              temperature; organ)
                        # BOTH accumulators, or the filling one could drift
                        # between the two paths unnoticed - and terminal heat
                        # reads only the second.
                        @test standalone === kernel
                    end
                end
            end
        end
    end
end

@testset "Leaf temperature raises exposure above air temperature" begin
    # The reason organ temperature stays selectable in this configuration: at
    # the gate cells air temperature loses 50-84% of the exposure hours, so the
    # two are different measurements of the same event, not a refinement.
    config = DiurnalConfig(; steps = 24,
                             shape = Agrocosm.diurnal_shape_code(:sinusoid))
    cft = _exposure_cft(Agrocosm.cft1; threshold = 28.0)
    air = _standalone(cft, config, 12.0, 14.0, 24.0)
    leaf = _standalone(cft, config, 12.0, 14.0, 24.0; organ = _organ(3.0))
    @test leaf > air
    # A canopy that intercepts nothing has no departure, so the energy balance
    # must fall away exactly rather than approximately.
    @test _standalone(cft, config, 12.0, 14.0, 24.0; organ = _organ(0.0)) === air
end

@testset "Exposure boundaries are exact, not merely small" begin
    config = DiurnalConfig(; steps = 24,
                             shape = Agrocosm.diurnal_shape_code(:flat))
    cft = _exposure_cft(Agrocosm.cft1; threshold = 20.0)

    # No day, no exposure. A negative daylength would otherwise give a negative
    # interval and a negative duration.
    @test _standalone(cft, config, 12.0, 0.0, 40.0) == 0
    @test _standalone(cft, config, 12.0, -1.0, 40.0) == 0

    # Flat shape with any range is the daily mean at every sub-step, so the
    # integral is the daylength times the smoothed exceedance - a closed form
    # this can be checked against.
    width = _T(Agrocosm.STERILITY_SMOOTHING_WIDTH)
    for temperature in (10.0, 19.5, 20.0, 21.0, 35.0)
        expected = _T(14) * Agrocosm.smooth_exceedance(_T(temperature) - _T(20), width)
        @test _standalone(cft, config, 12.0, 14.0, temperature) ≈ expected rtol = 1e-5
    end

    # Monotone in the daily mean: duration above a threshold cannot fall as the
    # whole temperature course is shifted up.
    sinusoid = DiurnalConfig(; steps = 24,
                               shape = Agrocosm.diurnal_shape_code(:sinusoid))
    series = [_standalone(cft, sinusoid, 12.0, 14.0, t) for t in 10.0:2.0:40.0]
    @test issorted(series)
    @test last(series) > first(series)
end

@testset "No forcing is a no-op, not a zeroing" begin
    # `nothing` has to leave the field alone rather than clear it: the call site
    # is unconditional, and in a run with sub-daily photosynthesis the
    # assimilation kernel is the writer. Clearing here would silently disable
    # the sink on exactly the configuration it was built for.
    crop = Agrocosm.init_crop(_T, 1, identity)
    state = test_model_state(crop)
    Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours .= _T(7.5)
    Agrocosm.heat_exposure!(
        _exposure_cft(Agrocosm.cft1), state, _T[14.0], _T[30.0], nothing,
    )
    @test Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours[1] == _T(7.5)
end

@testset "Both precisions agree on the same forcing" begin
    # The sink's rate will be calibrated against these hours, so a precision
    # dependence in them would be a precision dependence in the calibration.
    config32 = DiurnalConfig(; steps = 24,
                               shape = Agrocosm.diurnal_shape_code(:sinusoid))
    cft32 = _exposure_cft(Agrocosm.cft1; threshold = 24.0)
    single = _standalone(cft32, config32, 14.0, 14.0, 26.0; organ = _organ(3.0))

    base64 = convert_precision(Float64, Agrocosm.cft1)
    cft64 = Agrocosm.CFTParameters{Float64, Int32}(;
        (f => (f === :sterility_temperature ? 24.0 : getfield(base64, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
    crop = Agrocosm.init_crop(Float64, 1, identity)
    state = test_model_state(crop)
    organ = Agrocosm.OrganTemperatureForcing(
        [0.006], [101325.0], [2.0], [250.0], [-60.0], [0.2], [3.0], [3.0],
    )
    Agrocosm.heat_exposure!(
        cft64, state, [14.0], [26.0], DiurnalForcing(config32, [14.0]); organ,
    )
    @test Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours[1] ≈ single rtol = 1e-4
end

@testset "The closed form is the hard-threshold limit of the sub-daily integral" begin
    # The cell's claim is that it is the SAME integral taken analytically, so
    # away from the threshold it must converge to the numerical one as the
    # logistic separating them narrows.
    #
    # "Away from the threshold" is not a convenience: on a day whose maximum
    # sits AT the threshold the two can never agree at any width, because
    # `smooth_exceedance(0, w) = 0.5` for every w, so the entire day sits in the
    # transition. That degenerate case is asserted separately below, as the
    # finding it is rather than as a tolerance.
    clear(mean, range) = abs((mean + range / 2) - 35.0) > 2.0
    errors = Float64[]
    for width in (1.0, 0.3, 0.1, 0.03)
        worst = 0.0
        for mean in 20.0:2.0:40.0, range in 2.0:2.0:20.0
            clear(mean, range) || continue
            closed = Agrocosm.daily_statistic_exposure_hours(mean, range, 14.0, 35.0)
            numeric = Agrocosm.diurnal_heat_exposure(
                96, mean, range, 14.0, Agrocosm.DIURNAL_SINUSOID, 35.0, width)
            worst = max(worst, abs(closed - numeric))
        end
        push!(errors, worst)
    end
    @test issorted(errors; rev = true)
    @test last(errors) < 0.35
    # 1.38 h at the production width against 0.24 h at 0.03, measured on the
    # cases that are NOT near the threshold - so the smoothing is doing real
    # work even where the day is unambiguously hot or unambiguously not.
    @test first(errors) > 1.0

    # The smoothing finding, stated as an assertion so it cannot quietly change:
    # the production width credits a full 5 h of exposure to a day that never
    # crosses the threshold.
    peak_at_threshold = Agrocosm.diurnal_heat_exposure(
        96, 34.0, 2.0, 14.0, Agrocosm.DIURNAL_SINUSOID, 35.0, 1.0)
    @test peak_at_threshold > 5.0
    @test Agrocosm.daily_statistic_exposure_hours(34.0, 2.0, 14.0, 35.0) == 0

    # Boundaries, exactly.
    @test Agrocosm.daily_statistic_exposure_hours(40.0, 10.0, 0.0, 35.0) == 0
    @test Agrocosm.daily_statistic_exposure_hours(10.0, 4.0, 14.0, 35.0) == 0
    @test Agrocosm.daily_statistic_exposure_hours(50.0, 4.0, 14.0, 35.0) == 14

    # Monotone in the daily mean AND in the range, with no discontinuity as the
    # range closes - the defect that a partially smoothed threshold introduced.
    means = [Agrocosm.daily_statistic_exposure_hours(m, 10.0, 14.0, 35.0)
             for m in 25.0:1.0:45.0]
    ranges = [Agrocosm.daily_statistic_exposure_hours(32.0, r, 14.0, 35.0)
              for r in 0.0:0.5:20.0]
    hot_ranges = [Agrocosm.daily_statistic_exposure_hours(38.0, r, 14.0, 35.0)
                  for r in 0.0:0.5:20.0]
    @test issorted(means)
    @test issorted(ranges)
    # Above the threshold a wider range SHRINKS the duration, because the course
    # dips further below at night - so monotonicity has to reverse, and a flat
    # hot day must be the whole window.
    @test issorted(hot_ranges; rev = true)
    @test first(hot_ranges) == 14

    # The peak of the reconstructed course is the daily maximum, so a day whose
    # maximum sits exactly at the threshold has zero duration above it.
    @test Agrocosm.daily_statistic_exposure_hours(30.0, 10.0, 14.0, 35.0) == 0
    @test Agrocosm.daily_statistic_exposure_hours(30.0, 10.2, 14.0, 35.0) > 0
end

@testset "The daily-statistic kernel sees a mean-preserving range change" begin
    # THE point of the cell. Perturbation B holds the daily mean and widens the
    # range, and a criterion built from tasmax/tasmin is not blind to it - which
    # is why the paper cannot claim that every daily model responds zero.
    cft = _exposure_cft(Agrocosm.cft1; threshold = 30.0)
    function statistic(range)
        crop = Agrocosm.init_crop(_T, 1, identity)
        state = test_model_state(crop)
        Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours .= _T(-1)
        Agrocosm.daily_statistic_exposure!(
            cft, state, _T[14.0], _T[28.0], _T[range], true,
        )
        return Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours[1]
    end
    base = statistic(8.0)
    widened = statistic(20.0)
    @test widened > base
    @test base >= 0

    # And disabled is a no-op, not a zeroing, for the same reason as the
    # sub-daily pass: another kernel may be the writer.
    crop = Agrocosm.init_crop(_T, 1, identity)
    state = test_model_state(crop)
    Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours .= _T(3.25)
    Agrocosm.daily_statistic_exposure!(
        cft, state, _T[14.0], _T[28.0], _T[8.0], false,
    )
    @test Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours[1] == _T(3.25)
end

@testset "Both exposure kernels are per-cell independent" begin
    # Every other test in this file runs one cell, so a kernel indexing `[1]`
    # where it means `[cell]` would pass all of them and then produce one cell's
    # answer everywhere on a real domain - and silently, because cell 1 would be
    # right. This is the test that catches that class of bug, and it is the
    # property GPU parallelism needs: no cell may read another's inputs.
    cells = 5
    cft = _exposure_cft(Agrocosm.cft1; threshold = 26.0)
    config = DiurnalConfig(; steps = 24,
                             shape = Agrocosm.diurnal_shape_code(:sinusoid))
    temperatures = _T[18.0, 24.0, 28.0, 32.0, 38.0]
    ranges = _T[4.0, 8.0, 12.0, 16.0, 20.0]
    daylengths = _T[10.0, 12.0, 14.0, 16.0, 13.0]
    lais = _T[0.0, 1.0, 2.0, 3.0, 4.0]

    organ = Agrocosm.OrganTemperatureForcing(
        fill(_T(0.006), cells), fill(_T(101325.0), cells), fill(_T(2.0), cells),
        _T[150, 200, 250, 300, 350], fill(_T(-60.0), cells),
        fill(_T(0.2), cells), lais, fill(_T(3.0), cells),
    )

    crop = Agrocosm.init_crop(_T, cells, identity)
    state = test_model_state(crop)
    stress = Agrocosm.crop_stress_auxiliary(state)
    Agrocosm.heat_exposure!(
        cft, state, daylengths, temperatures,
        DiurnalForcing(config, ranges); organ,
    )
    batch_sterility = copy(stress.heat_exposure_hours)
    batch_filling = copy(stress.filling_exposure_hours)

    Agrocosm.daily_statistic_exposure!(
        cft, state, daylengths, temperatures, ranges, true,
    )
    batch_closed = copy(stress.heat_exposure_hours)

    # Each cell must equal what a one-cell run with that cell's inputs gives -
    # bitwise, because it is literally the same arithmetic.
    for index in 1:cells
        single_crop = Agrocosm.init_crop(_T, 1, identity)
        single = test_model_state(single_crop)
        single_organ = Agrocosm.OrganTemperatureForcing(
            _T[0.006], _T[101325.0], _T[2.0], _T[organ.shortwave[index]],
            _T[-60.0], _T[0.2], _T[lais[index]], _T[3.0],
        )
        Agrocosm.heat_exposure!(
            cft, single, _T[daylengths[index]], _T[temperatures[index]],
            DiurnalForcing(config, _T[ranges[index]]); organ = single_organ,
        )
        single_stress = Agrocosm.crop_stress_auxiliary(single)
        @test single_stress.heat_exposure_hours[1] === batch_sterility[index]
        @test single_stress.filling_exposure_hours[1] === batch_filling[index]

        Agrocosm.daily_statistic_exposure!(
            cft, single, _T[daylengths[index]], _T[temperatures[index]],
            _T[ranges[index]], true,
        )
        @test single_stress.heat_exposure_hours[1] === batch_closed[index]
    end

    # And the cells must actually differ, or the check above is vacuous.
    @test length(unique(batch_sterility)) == cells
    @test length(unique(batch_closed)) > 1
end
