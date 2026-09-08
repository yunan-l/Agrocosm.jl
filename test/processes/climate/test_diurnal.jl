using Agrocosm
using Test

# Gates G1-G4 of docs/01_subdaily_design.md. The scalar profile tests are
# self-contained; the kernel tests use the ordinary crop fixture.

const DIURNAL_STEPS = (1, 2, 3, 4, 6, 8, 12, 24, 48, 96)
const DIURNAL_DAYLENGTHS = (8.0, 10.0, 12.0, 14.5, 16.0)
const DIURNAL_RANGES = (0.0, 5.0, 12.0, 20.0, 25.0)
const DIURNAL_SHAPES = (DIURNAL_FLAT, DIURNAL_SINUSOID, DIURNAL_DAYTIME_NEUTRAL)

@testset "Diurnal radiation weights conserve the daily total (G3)" begin
    for T in (Float32, Float64)
        for shape in DIURNAL_SHAPES, steps in DIURNAL_STEPS
            total = sum(Agrocosm.diurnal_radiation_fraction(index, steps, shape, T)
                        for index in 1:steps)
            @test isapprox(total, one(T); atol = 32 * eps(T))
        end
        # At one step the fraction is exactly one, which is what makes the
        # sub-daily kernel reproduce `apar / daylength` bit for bit.
        for shape in DIURNAL_SHAPES
            @test Agrocosm.diurnal_radiation_fraction(1, 1, shape, T) === one(T)
        end
    end
end

@testset "Diurnal temperature conserves the daily mean (G4)" begin
    T = Float64
    mean_temperature = T(20)
    # :daytime_neutral - the sub-step mean is the daily mean exactly, any N.
    for daylength in DIURNAL_DAYLENGTHS, steps in DIURNAL_STEPS, range in DIURNAL_RANGES
        values = [Agrocosm.diurnal_temperature(index, steps, mean_temperature,
                                               T(range), T(daylength),
                                               DIURNAL_DAYTIME_NEUTRAL)
                  for index in 1:steps]
        @test isapprox(sum(values) / steps, mean_temperature; atol = 1e-11)
    end
    # :sinusoid - the continuous 24 h mean is the daily mean.
    for range in DIURNAL_RANGES
        samples = 200_000
        total = sum(mean_temperature +
                    T(range) * T(0.5) * sin(2 * pi * ((k + 0.5) * 24 / samples - 8) / 24)
                    for k in 0:(samples - 1))
        @test isapprox(total / samples, mean_temperature; atol = 1e-8)
    end
    # :flat ignores the range entirely.
    for daylength in DIURNAL_DAYLENGTHS, steps in DIURNAL_STEPS, range in DIURNAL_RANGES
        for index in 1:steps
            @test Agrocosm.diurnal_temperature(index, steps, mean_temperature,
                                               T(range), T(daylength),
                                               DIURNAL_FLAT) == mean_temperature
        end
    end
    # A zero range collapses every shape onto the daily mean.
    for shape in DIURNAL_SHAPES, daylength in DIURNAL_DAYLENGTHS, steps in DIURNAL_STEPS
        for index in 1:steps
            @test Agrocosm.diurnal_temperature(index, steps, mean_temperature,
                                               zero(T), T(daylength), shape) ==
                mean_temperature
        end
    end
end

@testset "Closed-form sub-step mean matches brute force" begin
    T = Float64
    for daylength in DIURNAL_DAYLENGTHS, steps in DIURNAL_STEPS
        brute = sum(Agrocosm.diurnal_shape_value(index, steps, T(daylength))
                    for index in 1:steps) / steps
        @test isapprox(brute, Agrocosm.diurnal_shape_substep_mean(steps, T(daylength));
                       atol = 1e-12)
    end
    # The single sub-step sits at solar noon.
    for daylength in DIURNAL_DAYLENGTHS
        @test Agrocosm.diurnal_substep_time(1, 1, T(daylength)) == T(12)
    end
end

@testset "Smooth exceedance is a usable threshold replacement" begin
    for T in (Float32, Float64)
        @test Agrocosm.smooth_exceedance(zero(T), T(0.5)) ≈ T(0.5)
        @test Agrocosm.smooth_exceedance(T(20), T(0.5)) ≈ one(T) atol = 1e-6
        @test Agrocosm.smooth_exceedance(T(-20), T(0.5)) ≈ zero(T) atol = 1e-6
        # Both tails must stay finite; a naive logistic overflows on one of them.
        @test isfinite(Agrocosm.smooth_exceedance(T(1e3), T(0.5)))
        @test isfinite(Agrocosm.smooth_exceedance(T(-1e3), T(0.5)))
        # A zero width recovers the hard indicator.
        @test Agrocosm.smooth_exceedance(T(0.1), zero(T)) == one(T)
        @test Agrocosm.smooth_exceedance(T(-0.1), zero(T)) == zero(T)
    end
end

@testset "Heat exposure is bounded and monotone" begin
    T = Float64
    daylength, range, critical, width = T(14), T(12), T(32), T(0.5)
    previous = -one(T)
    for mean_temperature in T.(20:2:40)
        hours = Agrocosm.diurnal_heat_exposure(48, mean_temperature, range, daylength,
                                               DIURNAL_SINUSOID, critical, width)
        @test zero(T) <= hours <= daylength
        @test hours > previous
        previous = hours
    end
    # No daylight means no exposure.
    @test Agrocosm.diurnal_heat_exposure(48, T(40), range, zero(T),
                                         DIURNAL_SINUSOID, critical, width) == zero(T)
    # A colder day than the threshold, with no diurnal range, gives none.
    @test Agrocosm.diurnal_heat_exposure(48, T(10), zero(T), daylength,
                                         DIURNAL_SINUSOID, critical, width) ≈ zero(T) atol = 1e-9
end

# --- Kernel-level gates -------------------------------------------------------

function _daily_stress(cft, daylength::T, temperature::T) where {T}
    return Agrocosm.compute_photosynthesis_temperature_stress(
        daylength, temperature, cft.path, cft.temp_co2, cft.temp_photos,
        T(Agrocosm.photoparams.tmc3), T(Agrocosm.photoparams.tmc4),
    )
end

_carbon_snapshot(crop) = (
    gross = copy(crop.fluxes.carbon.gross_assimilation),
    net = copy(crop.fluxes.carbon.net_assimilation),
    leaf = copy(crop.fluxes.carbon.leaf_respiration),
    water_limited = copy(crop.fluxes.carbon.water_limited_assimilation),
    vcmax = copy(crop.auxiliary.photosynthesis.vcmax),
)

@testset "Sub-daily C3 degenerates onto the daily kernel (G1, G2)" begin
    T = Float32
    apar, daylength, temperature, co2 = T[10.0], T[14.0], T[24.0], T[40.0]
    stress = _daily_stress(cft1, daylength[1], temperature[1])
    @test stress > T(1e-2)

    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.auxiliary.photosynthesis.temperature_stress .= stress
    photosynthesis_C3!(cft1, state, apar, daylength, temperature, co2;
                       comp_vcmax = true)
    daily = _carbon_snapshot(crop)
    # `_carbon_snapshot` keeps full arrays so the degeneracy checks below can
    # compare every field elementwise; index here for the scalar comparison.
    @test daily.gross[1] > zero(T)

    # One step with a zero range reproduces the daily result exactly, whatever
    # the shape. `nothing` must route to the daily kernel identically.
    for shape in DIURNAL_SHAPES
        crop = init_crop(1, identity)
        state = test_model_state(crop)
        crop.auxiliary.photosynthesis.temperature_stress .= stress
        photosynthesis!(Val(:C3), cft1, state, apar, daylength, temperature, co2,
                        DiurnalForcing(T[0.0]; steps = 1, shape = shape);
                        comp_vcmax = true)
        current = _carbon_snapshot(crop)
        for field in keys(daily)
            @test getproperty(current, field) == getproperty(daily, field)
        end
    end

    # One step also reproduces it exactly at a non-zero range for :flat and
    # :daytime_neutral, because their single sub-step sits at the daily mean.
    for shape in (DIURNAL_FLAT, DIURNAL_DAYTIME_NEUTRAL)
        crop = init_crop(1, identity)
        state = test_model_state(crop)
        crop.auxiliary.photosynthesis.temperature_stress .= stress
        photosynthesis!(Val(:C3), cft1, state, apar, daylength, temperature, co2,
                        DiurnalForcing(T[14.0]; steps = 1, shape = shape);
                        comp_vcmax = true)
        @test _carbon_snapshot(crop).gross == daily.gross
    end

    # Passing `nothing` is the daily path.
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.auxiliary.photosynthesis.temperature_stress .= stress
    photosynthesis!(Val(:C3), cft1, state, apar, daylength, temperature, co2,
                    nothing; comp_vcmax = true)
    @test _carbon_snapshot(crop).gross == daily.gross
end

@testset "Sub-daily integration lowers daily assimilation (Jensen sign)" begin
    # Float32 throughout: `init_crop(1, identity)` builds Float32 state and the
    # kernels require one element type across every array argument.
    T = Float32
    apar, daylength, co2 = T[10.0], T[14.0], T[40.0]

    function subdaily_gross(cft, pathway, temperature, range, steps, shape)
        crop = init_crop(1, identity)
        state = test_model_state(crop)
        crop.auxiliary.photosynthesis.temperature_stress .=
            _daily_stress(cft, daylength[1], temperature)
        photosynthesis!(pathway, cft, state, apar, daylength, T[temperature], co2,
                        DiurnalForcing(T[range]; steps = steps, shape = shape);
                        comp_vcmax = true)
        return crop.fluxes.carbon.gross_assimilation[1]
    end

    for (cft, pathway) in ((cft1, Val(:C3)), (cft3, Val(:C4)))
        reference = subdaily_gross(cft, pathway, T(24), T(0), 1, DIURNAL_SINUSOID)
        @test reference > zero(T)

        # Channel 1: radiation curvature alone. With a zero diurnal range the
        # temperature never moves, yet redistributing the day's PAR onto a
        # half-sine cannot raise the total, because co-limited assimilation is
        # concave in the light-limited rate and the sub-step mean of that rate
        # is unchanged (Jensen).
        #
        # The inequality is weak on purpose. `theta = 0.99` makes co-limitation
        # nearly `min(je, jc)`, i.e. almost piecewise linear, so the effect is
        # only strict when the sub-steps straddle the je = jc corner. A day that
        # is firmly light- or Rubisco-limited has no radiation curvature at all.
        # Asserting the guaranteed direction avoids a test that depends on where
        # this particular fixture happens to sit.
        radiation_only = subdaily_gross(cft, pathway, T(24), T(0), 48, DIURNAL_SINUSOID)
        @test radiation_only <= reference * (one(T) + T(1e-6))

        # Channel 2: temperature curvature, at a hot daily mean, with the
        # sub-step mean temperature pinned to the daily mean by
        # :daytime_neutral. This isolates curvature from daytime warming.
        hot_reference = subdaily_gross(cft, pathway, T(34), T(0), 48, DIURNAL_SINUSOID)
        hot_curvature = subdaily_gross(cft, pathway, T(34), T(12), 48,
                                       DIURNAL_DAYTIME_NEUTRAL)
        @test hot_curvature < hot_reference

        # Channel 3: curvature plus the real daytime warming. Measured to be
        # the larger of the two for wheat and maize.
        hot_physical = subdaily_gross(cft, pathway, T(34), T(12), 48, DIURNAL_SINUSOID)
        @test hot_physical < hot_curvature

        # Refining the step count must converge rather than drift. Compared as
        # a relative change so Float32 round-off at small differences cannot
        # decide the test.
        coarse = subdaily_gross(cft, pathway, T(34), T(12), 24, DIURNAL_SINUSOID)
        fine = subdaily_gross(cft, pathway, T(34), T(12), 96, DIURNAL_SINUSOID)
        finer = subdaily_gross(cft, pathway, T(34), T(12), 192, DIURNAL_SINUSOID)
        scale = max(abs(hot_physical), eps(T))
        @test abs(finer - fine) / scale <= abs(fine - coarse) / scale + T(1e-5)
        @test abs(finer - fine) / scale < T(0.02)
    end
end
