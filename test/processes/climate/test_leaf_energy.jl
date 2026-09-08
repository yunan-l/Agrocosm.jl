using Agrocosm
using Agrocosm.KernelAbstractions   # reached through Agrocosm, not a test dependency
using Test
using JLD2

# Gates E1-E6 of docs/04_organ_temperature_design.md. E7 (Enzyme) lives with the
# rest of the AD suite in test/ad/test_weather_attribution.jl.
#
# The forcing is the real single-cell wheat series shipped in
# examples/climate_2000_2009.jld2 (10 years of daily temperature, shortwave and
# net longwave). That archive predates this module and carries no wind, specific
# humidity or surface pressure, so those three use documented placeholder
# constants below; the server's GSWP3-W5E5 supplies the real fields
# (`sfcwind`, `huss`, `ps`) once the runner side is wired.

const EXAMPLE_CLIMATE_PATH =
    normpath(joinpath(@__DIR__, "..", "..", "..", "examples", "climate_2000_2009.jld2"))

# Placeholders, not measurements. Chosen to be unremarkable for a temperate
# rainfed cell so that any extreme behaviour in a test comes from the swept
# variable rather than from these.
const PLACEHOLDER_WIND = 2.0          # m s-1
const PLACEHOLDER_PRESSURE = 101325.0 # Pa, sea level
const PLACEHOLDER_HUMIDITY = 0.006    # kg kg-1, ~70% RH at 10 C

const LEAF_DIMENSION = 0.04           # m, CLM5's uniform value
const LEAF_EMISSIVITY = 0.98
const EXTINCTION = 0.5                # cft1/cft3 lightextcoeff

"""Solve the full non-linear balance by bisection, for use as ground truth."""
function exact_leaf_temperature(air, shortwave, longwave, humidity, pressure,
                                wind, conductance, lai;
                                dimension = LEAF_DIMENSION,
                                emissivity = LEAF_EMISSIVITY,
                                extinction = EXTINCTION)
    residual(t) = Agrocosm.leaf_energy_residual(
        t, air, shortwave, longwave, humidity, pressure, wind, conductance, lai,
        dimension, emissivity, extinction,
    )
    low, high = air - 40.0, air + 40.0
    residual(low) * residual(high) > 0 && return NaN
    for _ in 1:200
        mid = (low + high) / 2
        residual(low) * residual(mid) <= 0 ? (high = mid) : (low = mid)
    end
    return (low + high) / 2
end

departure(air, shortwave, longwave, humidity, pressure, wind, conductance, lai;
          dimension = LEAF_DIMENSION, emissivity = LEAF_EMISSIVITY,
          extinction = EXTINCTION) = Agrocosm.leaf_temperature_departure(
    air, shortwave, longwave, humidity, pressure, wind, conductance, lai,
    dimension, emissivity, extinction,
)

@testset "Saturation curve matches the model's existing one" begin
    # radiation.jl's compute_equilibrium_evaporation embeds the same Tetens
    # coefficients; its `slope` term must equal our analytic derivative, or the
    # model would carry two different saturation curves.
    for temperature in (-5.0, 0.0, 10.0, 25.0, 35.0, 45.0)
        offset = 237.3 + temperature
        embedded = 2.503e6 * exp(17.269 * temperature / offset) / offset^2
        @test Agrocosm.saturation_vapour_pressure_slope(temperature) ≈ embedded rtol = 1e-4
        # And the derivative really is the derivative of the curve.
        h = 1e-6
        numeric = (Agrocosm.saturation_vapour_pressure(temperature + h) -
                   Agrocosm.saturation_vapour_pressure(temperature - h)) / (2h)
        @test Agrocosm.saturation_vapour_pressure_slope(temperature) ≈ numeric rtol = 1e-6
    end
    @test Agrocosm.saturation_vapour_pressure(0.0) ≈ 610.78
end

@testset "Vapour pressure from specific humidity" begin
    # Round trip against the textbook inverse q = 0.622 e / (p - 0.378 e).
    for q in (0.001, 0.006, 0.015, 0.03), p in (70000.0, 101325.0)
        e = Agrocosm.vapour_pressure_from_specific_humidity(q, p)
        @test 0.622 * e / (p - 0.378 * e) ≈ q rtol = 1e-10
    end
    @test Agrocosm.vapour_pressure_from_specific_humidity(0.0, 101325.0) == 0.0
end

@testset "E2: departure vanishes as the boundary layer thins" begin
    # A thinner boundary layer welds the leaf to the air. The sweep stops at
    # 1e-3 m because `leaf_boundary_layer_conductance` floors the characteristic
    # dimension there to keep the square root away from zero, so this limit is
    # asymptotic, not exact: at the floor a residual departure of order 0.1 K
    # survives. The ablation experiment's exact "organ temperature off" arm
    # therefore has to come from the configuration switch, never from shrinking
    # this parameter.
    previous = Inf
    for dimension in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3)
        value = abs(departure(30.0, 400.0, -60.0, PLACEHOLDER_HUMIDITY,
                              PLACEHOLDER_PRESSURE, PLACEHOLDER_WIND, 3.0, 3.0;
                              dimension))
        @test value < previous
        previous = value
    end
    @test previous < 0.2

    # Bare ground intercepts nothing. This degeneracy *is* exact.
    @test departure(30.0, 400.0, -60.0, PLACEHOLDER_HUMIDITY, PLACEHOLDER_PRESSURE,
                    PLACEHOLDER_WIND, 3.0, 0.0) == 0.0
end

@testset "E3: closed form satisfies the full non-linear balance" begin
    # Sweep the extreme corner: hot, dry, bright, low wind, near-closed stomata.
    worst_residual = 0.0
    worst_error = 0.0
    for air in (25.0, 35.0, 45.0), shortwave in (300.0, 500.0, 700.0),
        conductance in (0.5, 2.0, 8.0), wind in (0.5, 2.0, 4.0),
        humidity in (0.003, 0.008, 0.015), lai in (1.0, 3.0, 6.0)

        estimate = departure(air, shortwave, -60.0, humidity, PLACEHOLDER_PRESSURE,
                             wind, conductance, lai)
        @test isfinite(estimate)
        residual = Agrocosm.leaf_energy_residual(
            air + estimate, air, shortwave, -60.0, humidity, PLACEHOLDER_PRESSURE,
            wind, conductance, lai, LEAF_DIMENSION, LEAF_EMISSIVITY, EXTINCTION,
        )
        worst_residual = max(worst_residual, abs(residual))

        truth = exact_leaf_temperature(air, shortwave, -60.0, humidity,
                                       PLACEHOLDER_PRESSURE, wind, conductance, lai)
        isnan(truth) && continue
        worst_error = max(worst_error, abs(air + estimate - truth))
    end
    @info "E3 closed-form accuracy" worst_residual worst_error
    # The temperature threshold is the meaningful one: 0.05 K is an order of
    # magnitude below the scheme's own structural biases (no stability
    # correction, conductance frozen across sub-steps), so tightening it further
    # would be measuring noise against a much larger systematic error. The flux
    # threshold is the same statement rescaled - a 0.05 K miss against a typical
    # denominator of ~60 W m-2 K-1 is ~3 W m-2, itself under 1% of the several
    # hundred W m-2 of net radiation being balanced.
    @test worst_error < 0.05       # K
    @test worst_residual < 3.0     # W m-2
end

@testset "E4: the re-expansion earns its place" begin
    # One expansion about air temperature vs. the shipped two-pass version,
    # both against the bisection solution. The design claims two orders of
    # magnitude; assert one and report the rest.
    worst_single = 0.0
    worst_double = 0.0
    for air in (30.0, 38.0, 45.0), shortwave in (500.0, 700.0),
        conductance in (0.3, 1.0), wind in (0.5, 1.5), humidity in (0.003, 0.006)

        truth = exact_leaf_temperature(air, shortwave, -60.0, humidity,
                                       PLACEHOLDER_PRESSURE, wind, conductance, 3.0)
        isnan(truth) && continue

        vapour = Agrocosm.vapour_pressure_from_specific_humidity(humidity, PLACEHOLDER_PRESSURE)
        density = Agrocosm.molar_air_density(air, PLACEHOLDER_PRESSURE)
        area = 2 * Agrocosm.effective_leaf_area(3.0, EXTINCTION)
        heat = area * Agrocosm.leaf_boundary_layer_conductance(wind, LEAF_DIMENSION)
        vapour_conductance = Agrocosm.series_conductance(
            conductance * 1e-3 * density, heat * (0.147 / 0.135),
        )
        cover = Agrocosm.canopy_cover_fraction(3.0, EXTINCTION)
        single = Agrocosm._departure_about(
            air, air, shortwave, -60.0, vapour, PLACEHOLDER_PRESSURE,
            vapour_conductance, heat, LEAF_EMISSIVITY, cover,
        )
        double = departure(air, shortwave, -60.0, humidity, PLACEHOLDER_PRESSURE,
                           wind, conductance, 3.0)

        worst_single = max(worst_single, abs(air + single - truth))
        worst_double = max(worst_double, abs(air + double - truth))
    end
    @info "E4 re-expansion" worst_single worst_double ratio = worst_single / worst_double
    @test worst_double < worst_single / 10
    @test worst_double < 0.01
end

@testset "E5: signs and monotonicity" begin
    base = (35.0, 500.0, -60.0, PLACEHOLDER_PRESSURE, 3.0)
    air, shortwave, longwave, pressure, lai = base

    # Moister air means a smaller deficit, so less evaporative cooling is
    # available and the leaf sits warmer.
    humid = [departure(air, shortwave, longwave, q, pressure, PLACEHOLDER_WIND, 3.0, lai)
             for q in (0.002, 0.006, 0.012, 0.020)]
    @test issorted(humid)

    # Closing stomata removes transpirational cooling: hotter.
    closing = [departure(air, shortwave, longwave, PLACEHOLDER_HUMIDITY, pressure,
                         PLACEHOLDER_WIND, g, lai) for g in (12.0, 6.0, 2.0, 0.5)]
    @test issorted(closing)

    # More wind couples the leaf to the air, whatever the sign of the departure.
    for conductance in (0.5, 10.0)
        magnitudes = [abs(departure(air, shortwave, longwave, PLACEHOLDER_HUMIDITY,
                                    pressure, u, conductance, lai))
                      for u in (0.5, 1.0, 2.0, 4.0, 8.0)]
        @test issorted(magnitudes; rev = true)
    end

    # More radiation loaded onto the canopy: hotter.
    brighter = [departure(air, s, longwave, PLACEHOLDER_HUMIDITY, pressure,
                          PLACEHOLDER_WIND, 2.0, lai) for s in (100.0, 300.0, 500.0, 800.0)]
    @test issorted(brighter)
end

@testset "E6: the irrigation sign flip is reproduced" begin
    # The observational signature this scheme has to reproduce: irrigated rye
    # measured near -2 C, rainfed rye on sandy soil up to +7.5 C
    # (Siebert et al. 2014, ERL 9:044012). Same equation, both signs, driven
    # only by conductance.
    hot_dry_day = (35.0, 600.0, -60.0, 0.006, PLACEHOLDER_PRESSURE, 1.5, 3.0)
    air, shortwave, longwave, humidity, pressure, wind, lai = hot_dry_day

    well_watered = departure(air, shortwave, longwave, humidity, pressure, wind, 15.0, lai)
    water_limited = departure(air, shortwave, longwave, humidity, pressure, wind, 0.5, lai)

    @info "E6 irrigation sign flip" well_watered water_limited
    @test well_watered < 0        # transpirational cooling
    @test water_limited > 0       # canopy runs hot
    @test water_limited - well_watered > 3.0
end

@testset "Real single-cell forcing from examples/" begin
    @test isfile(EXAMPLE_CLIMATE_PATH)
    climate = load(EXAMPLE_CLIMATE_PATH)["climate"]
    # Column 1 is the first of the ten example cells; temp/swdown/lwnet are real
    # daily series, the other three inputs are the placeholders declared above.
    temperature = Float64.(climate.temp[:, 1])
    shortwave = Float64.(climate.swdown[:, 1])
    longwave = Float64.(climate.lwnet[:, 1])
    @test length(temperature) == 3650

    for conductance in (0.5, 3.0, 12.0)
        departures = [
            departure(temperature[day], shortwave[day], longwave[day],
                      PLACEHOLDER_HUMIDITY, PLACEHOLDER_PRESSURE, PLACEHOLDER_WIND,
                      conductance, 3.0)
            for day in eachindex(temperature)
        ]
        @test all(isfinite, departures)
        # Nothing physical on a temperate cell should leave this envelope; a
        # unit error of the kind that multiplies conductance by molar density
        # would either flatten these to zero or blow past it.
        @test all(d -> -15.0 < d < 25.0, departures)

        # Every day must still satisfy the balance it was solved from.
        worst = maximum(
            abs(Agrocosm.leaf_energy_residual(
                temperature[day] + departures[day], temperature[day], shortwave[day],
                longwave[day], PLACEHOLDER_HUMIDITY, PLACEHOLDER_PRESSURE,
                PLACEHOLDER_WIND, conductance, 3.0, LEAF_DIMENSION,
                LEAF_EMISSIVITY, EXTINCTION,
            )) for day in eachindex(temperature)
        )
        @info "Real-forcing residual" conductance worst extremes = extrema(departures)
        @test worst < 1.0
    end
end

@kernel inbounds = true function leaf_departure_kernel!(
    output::AbstractVector{T},
    air::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    humidity::AbstractVector{T},
    pressure::AbstractVector{T},
    wind::AbstractVector{T},
    conductance::AbstractVector{T},
    lai::AbstractVector{T},
    dimension::T,
    emissivity::T,
    extinction::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    output[cell] = Agrocosm.leaf_temperature_departure(
        air[cell], shortwave[cell], longwave[cell], humidity[cell], pressure[cell],
        wind[cell], conductance[cell], lai[cell], dimension, emissivity, extinction,
    )
end

@testset "Runs inside a KernelAbstractions kernel" begin
    # The scalar functions are only useful if they survive a kernel launch: no
    # allocation, no dynamic dispatch, no Float64 sneaking in through an
    # unconverted global constant. Launching on the CPU backend exercises the
    # same code path the GPU backend compiles, which is what makes the
    # single-precision result below meaningful.
    climate = load(EXAMPLE_CLIMATE_PATH)["climate"]
    cells = 64
    T = Float32
    air = T.(climate.temp[1:cells, 1])
    shortwave = T.(climate.swdown[1:cells, 1])
    longwave = T.(climate.lwnet[1:cells, 1])
    humidity = fill(T(PLACEHOLDER_HUMIDITY), cells)
    pressure = fill(T(PLACEHOLDER_PRESSURE), cells)
    wind = fill(T(PLACEHOLDER_WIND), cells)
    conductance = fill(T(3), cells)
    lai = fill(T(3), cells)
    output = zeros(T, cells)

    Agrocosm.launch_1D!(
        leaf_departure_kernel!, output, air, shortwave, longwave, humidity,
        pressure, wind, conductance, lai, T(LEAF_DIMENSION), T(LEAF_EMISSIVITY),
        T(EXTINCTION),
    )

    @test eltype(output) === Float32
    @test all(isfinite, output)
    for cell in 1:cells
        scalar = Agrocosm.leaf_temperature_departure(
            air[cell], shortwave[cell], longwave[cell], humidity[cell],
            pressure[cell], wind[cell], conductance[cell], lai[cell],
            T(LEAF_DIMENSION), T(LEAF_EMISSIVITY), T(EXTINCTION),
        )
        @test output[cell] === scalar
    end
end

@testset "E1: organ temperature off is bitwise step 1" begin
    # The sub-daily kernel gained an argument and its loop body was rewritten,
    # so "logically equivalent" is not good enough: the switched-off path has to
    # reproduce step 1 bit for bit, and the degenerate switched-on path has to
    # match it too. Otherwise the ablation experiment's control arm is not a
    # control.
    T = Float32
    apar, daylength, temperature, co2 = T[10.0], T[14.0], T[30.0], T[40.0]
    range = T[12.0]
    config = DiurnalConfig(; steps = 24, shape = DIURNAL_SINUSOID)
    stress = T(0.8)

    function run(pathway, cft, organ)
        crop = init_crop(1, identity)
        state = test_model_state(crop)
        crop.auxiliary.photosynthesis.temperature_stress .= stress
        photosynthesis!(pathway, cft, state, apar, daylength, temperature, co2,
                        DiurnalForcing(config, range); comp_vcmax = true, organ)
        return copy(crop.fluxes.carbon.gross_assimilation)
    end

    # lai = 0: the canopy intercepts nothing, the departure is exactly zero, and
    # the whole energy balance has to fall away without leaving a trace.
    inert = Agrocosm.OrganTemperatureForcing(
        T[0.006], T[101325.0], T[2.0], T[250.0], T[-60.0], T[0.2], T[0.0], T[3.0],
    )
    # A real canopy, to prove the wiring is actually live.
    active = Agrocosm.OrganTemperatureForcing(
        T[0.006], T[101325.0], T[2.0], T[250.0], T[-60.0], T[0.2], T[3.0], T[3.0],
    )

    for (pathway, cft) in ((Val(:C3), cft1), (Val(:C4), cft3))
        off = run(pathway, cft, nothing)
        @test run(pathway, cft, inert) == off
        @test run(pathway, cft, active) != off
    end
end

@testset "Float32 agrees with Float64 on the real series" begin
    climate = load(EXAMPLE_CLIMATE_PATH)["climate"]
    temperature = climate.temp[:, 1]      # already Float32
    shortwave = climate.swdown[:, 1]
    longwave = climate.lwnet[:, 1]
    worst = 0.0f0
    for day in eachindex(temperature)
        single = Agrocosm.leaf_temperature_departure(
            temperature[day], shortwave[day], longwave[day],
            Float32(PLACEHOLDER_HUMIDITY), Float32(PLACEHOLDER_PRESSURE),
            Float32(PLACEHOLDER_WIND), 3.0f0, 3.0f0,
            Float32(LEAF_DIMENSION), Float32(LEAF_EMISSIVITY), Float32(EXTINCTION),
        )
        double = departure(Float64(temperature[day]), Float64(shortwave[day]),
                           Float64(longwave[day]), PLACEHOLDER_HUMIDITY,
                           PLACEHOLDER_PRESSURE, PLACEHOLDER_WIND, 3.0, 3.0)
        @test isfinite(single)
        worst = max(worst, abs(single - Float32(double)))
    end
    @info "Float32 vs Float64 on real forcing" worst
    @test worst < 1.0f-2
end
