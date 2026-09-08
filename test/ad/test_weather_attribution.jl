using Agrocosm
using Enzyme
using Test

include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

@testset "Prepared weather kernel" begin
    weather = Agrocosm.init_weather(Float64, 2, identity)
    forcing = reshape(collect(1.0:30.0), 3, 2, 5)
    Agrocosm.apply_weather_forcing!(weather, forcing, 2)
    for (variable, field) in enumerate((:temp, :prec, :swr, :lwr, :wind))
        @test getproperty(weather, field) == forcing[2, :, variable]
    end
    @test weather.no3_deposition == zeros(2)
    @test weather.annual_co2 == zeros(1)
end

for cft_id in (isempty(ARGS) ? (1, 3) : Tuple(parse.(Int, ARGS)))
    @testset "CFT $cft_id fixed-harvest weather AD" begin
        T = get(ENV, "WEATHER_TEST_PRECISION", "64") == "32" ? Float32 : Float64
        window_days = get(ENV, "WEATHER_TEST_WINDOW", "8") == "season" ? :season : 8
        sowing_day = parse(Int, get(ENV, "WEATHER_TEST_SOWING_DAY", "100"))
        case = weather_attribution_fixture(cft_id; T, window_days, sowing_day)
        if get(ENV, "WEATHER_TEST_DAILY_CO2", "0") == "1"
            # The server case builder uses the climate reader's daily CO₂
            # contract, not the annual-vector contract of older AD fixtures.
            daily = repeat(case.climate.co2; inner = 365)
            case = merge(case, (; climate = merge(case.climate, (co2 = daily, co2_daily = true))))
        end
        (; forcing, state, cft, parameters, climate, days, harvest_day) = case
        original = deepcopy(forcing)
        reference = weather_harvest_replay(forcing, state, cft, parameters, climate, days, harvest_day)
        @test reference.schedule_matches
        @test reference.yield > 0
        @info "Weather attribution production fixture" cft_id harvest_day production_yield = reference.yield
        flush(stderr)

        result = enzyme_weather_harvest_gradient(forcing, state, cft, parameters,
            climate, days, harvest_day; block_days = 4)
        @test result.primal ≈ result.production_yield rtol = 1e-3 atol = 1e-5
        @test result.reverse_primal ≈ result.primal rtol = 1e-10 atol = 1e-12
        @test all(isfinite, result.gradient)
        @test forcing == original
        @test all(iszero, result.gradient[1:(first(days) - 1), :, :])
        @test all(iszero, result.gradient[harvest_day:end, :, :])
        @test any(!iszero, result.gradient[days, :, 1])
        @test any(!iszero, result.gradient[days, :, 3])

        direction = zeros(T, size(forcing))
        direction[first(days):(first(days) + 2), 1, 1] .= 1.0
        direction[first(days):(first(days) + 2), 1, 2] .= 0.3
        projection = sum(result.gradient .* direction)

        # Float32 production includes a finite-precision bisection solve;
        # a 0.01 C perturbation can be below its useful differencing scale.
        # A 0.1/0.05/0.02/0.01 sweep validates 0.05 without relaxing tolerances.
        epsilon = T === Float32 ? T(0.05) : T(0.01)
        plus = weather_harvest_replay(forcing .+ epsilon .* direction, state,
            cft, parameters, climate, days, harvest_day)
        minus = weather_harvest_replay(forcing .- epsilon .* direction, state,
            cft, parameters, climate, days, harvest_day)
        @test plus.schedule_matches && minus.schedule_matches
        fd = (plus.yield - minus.yield) / (2epsilon)
        @info "Weather directional validation" cft_id ad = projection fd
        flush(stderr)
        @test projection ≈ fd rtol = 0.03 atol = 2e-5

        forward = enzyme_weather_forward_directional(forcing, direction, state,
            cft, parameters, climate, days, harvest_day)
        @test forward.directional ≈ projection rtol = 2e-5 atol = 1e-8

        full = enzyme_weather_harvest_gradient(forcing, state, cft, parameters,
            climate, days, harvest_day; block_days = length(days))
        @test full.gradient ≈ result.gradient rtol = 2e-5 atol = 1e-8

        # Separate channels avoid a temperature signal masking a missing
        # precipitation, radiation or wind derivative in the combined probe.
        for variable in 2:5
            probe = zeros(T, size(forcing))
            probe[first(days):(first(days) + 2), 1, variable] .= one(T)
            # Radiation is in W/m², not degrees or mm: use its independently
            # checked differencing scale (0.05 W/m² was noisy in Float32).
            delta = variable in (3, 4) ? T(0.1) : T(0.05)
            upper = weather_harvest_replay(forcing .+ delta .* probe, state,
                cft, parameters, climate, days, harvest_day)
            lower = weather_harvest_replay(forcing .- delta .* probe, state,
                cft, parameters, climate, days, harvest_day)
            @test upper.schedule_matches && lower.schedule_matches
            numeric = (upper.yield - lower.yield) / (2delta)
            @info "Weather channel validation" cft_id variable delta ad = sum(result.gradient .* probe) fd = numeric
            @test sum(result.gradient .* probe) ≈ numeric rtol = 0.03 atol = 2e-5
        end

        @test_throws ArgumentError enzyme_weather_harvest_gradient(forcing, state,
            cft, parameters, climate, days, harvest_day + 1)
        invalid = copy(forcing)
        invalid[first(days), 1, 2] = -1
        @test_throws ArgumentError weather_harvest_replay(invalid, state,
            cft, parameters, climate, days, harvest_day)
    end
end

for cft_id in (isempty(ARGS) ? (1, 3) : Tuple(parse.(Int, ARGS)))
    @testset "CFT $cft_id fixed-harvest weather AD with sub-daily photosynthesis" begin
        T = get(ENV, "WEATHER_TEST_PRECISION", "64") == "32" ? Float32 : Float64
        sowing_day = parse(Int, get(ENV, "WEATHER_TEST_SOWING_DAY", "100"))
        diurnal_config = DiurnalConfig(; steps = 24, shape = DIURNAL_SINUSOID)

        # G1/G2 in a real Enzyme pipeline: the sub-daily path must degenerate
        # exactly onto the daily kernel in every case the design says it does,
        # both in the ordinary replay and in the gradient that feeds
        # `enzyme_weather_harvest_gradient`. The exact-degeneracy rows are
        # tabulated in `src/processes/crop/photosynthesis_subdaily.jl`: a zero
        # diurnal range under any shape, and `:flat` or `:daytime_neutral`
        # under any range. `:sinusoid` with a nonzero range is deliberately
        # NOT one of them -- a single sub-step of a sinusoid sits at the
        # solar-noon temperature, not the daily mean -- so it is asserted
        # below as documented behaviour instead. This is checked end-to-end
        # through the AD path, which had no `diurnal_config` wiring at all
        # before 2026-09-08.
        daily_case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day)
        daily_reference = weather_harvest_replay(daily_case.forcing, daily_case.state,
            daily_case.cft, daily_case.parameters, daily_case.climate, daily_case.days,
            daily_case.harvest_day)
        subdaily_replay = (config; amplitude = T(10)) -> begin
            replay_case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day,
                diurnal_config = config, diurnal_amplitude = amplitude)
            replay = weather_harvest_replay(replay_case.forcing, replay_case.state,
                replay_case.cft, replay_case.parameters, replay_case.climate,
                replay_case.days, replay_case.harvest_day; diurnal_config = config)
            return (replay_case, replay)
        end

        degenerate_config = DiurnalConfig(; steps = 1, shape = DIURNAL_FLAT)
        degenerate_case, degenerate_reference = subdaily_replay(degenerate_config)
        @test degenerate_reference.yield == daily_reference.yield
        # Many sub-steps, flat shape: the radiation weights must sum to one and
        # every sub-step sees the daily mean. This is the row that actually
        # catches a mis-normalized integration weight, which a single-sub-step
        # check cannot see. Equality here is up to round-off, not bitwise:
        # summing `steps` contributions is not associative in floating point,
        # so the tolerance is a few hundred eps rather than zero. A real
        # normalization error is O(1), nowhere near this band.
        @test isapprox(last(subdaily_replay(DiurnalConfig(; steps = 24, shape = DIURNAL_FLAT))).yield,
            daily_reference.yield; rtol = 1000 * eps(T))
        # A zero diurnal range removes the temperature spread, so at one
        # sub-step every shape collapses onto the daily state exactly.
        @test last(subdaily_replay(DiurnalConfig(; steps = 1, shape = DIURNAL_SINUSOID);
            amplitude = T(0))).yield == daily_reference.yield
        # But a zero range does NOT make a non-flat shape degenerate once there
        # is more than one sub-step: the shape still redistributes the day's
        # PAR across sub-steps, and co-limited assimilation is concave in
        # light, so integrating it over an uneven light course gives less than
        # evaluating it once at the mean. That is the sub-daily scheme's second
        # Jensen channel -- light curvature -- and it is active independently
        # of the temperature curvature the diurnal range drives. Asserted
        # explicitly so the two channels cannot be silently conflated.
        light_only = last(subdaily_replay(DiurnalConfig(; steps = 24, shape = DIURNAL_SINUSOID);
            amplitude = T(0))).yield
        @test light_only != daily_reference.yield
        @test light_only < daily_reference.yield
        # `:daytime_neutral` subtracts its own closed-form sub-step mean, which
        # at one sub-step is the solar-noon value itself.
        @test last(subdaily_replay(DiurnalConfig(; steps = 1, shape = DIURNAL_DAYTIME_NEUTRAL))).yield ==
            daily_reference.yield
        # The documented non-degenerate row, pinned so a future change to the
        # shape machinery cannot silently turn it into degeneracy: one sub-step
        # of `:sinusoid` integrates at the solar-noon temperature, so it must
        # differ from the daily kernel while remaining a valid harvest.
        noon_reference = last(subdaily_replay(DiurnalConfig(; steps = 1, shape = DIURNAL_SINUSOID)))
        @test noon_reference.yield != daily_reference.yield
        @test isfinite(noon_reference.yield) && noon_reference.yield > zero(T)
        daily_gradient = enzyme_weather_harvest_gradient(daily_case.forcing, daily_case.state,
            daily_case.cft, daily_case.parameters, daily_case.climate, daily_case.days,
            daily_case.harvest_day; block_days = 4)
        degenerate_gradient = enzyme_weather_harvest_gradient(degenerate_case.forcing, degenerate_case.state,
            degenerate_case.cft, degenerate_case.parameters, degenerate_case.climate,
            degenerate_case.days, degenerate_case.harvest_day; block_days = 4, diurnal_config = degenerate_config)
        # The forward yield is bitwise identical (checked above), but the
        # gradient is not: Enzyme differentiates two different source-code
        # shapes here (the daily kernel's direct expression vs. the sub-daily
        # kernel's single-iteration loop), and the resulting adjoint code sums
        # floating-point contributions in a different order. That is a
        # non-associativity effect, not a modelling difference, so it is
        # bounded by a handful of ULP rather than zero -- confirmed directly:
        # of 1825 entries exactly one differs, by a relative 1.3e-7, which is
        # eps(Float32) itself. A real wiring bug would be orders of magnitude
        # larger.
        @test isapprox(degenerate_gradient.gradient, daily_gradient.gradient;
            rtol = 1000 * eps(T), atol = 1000 * eps(T))

        # The actual sub-daily case: 24 steps, a physically sized diurnal
        # range. This is the first time this repository differentiates
        # through `_enzyme_continuous_transition!`'s new `diurnal_config`
        # branch instead of just its forward replay.
        case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day,
            diurnal_config, diurnal_amplitude = T(10))
        (; forcing, state, cft, parameters, climate, days, harvest_day) = case
        original = deepcopy(forcing)
        reference = weather_harvest_replay(forcing, state, cft, parameters, climate,
            days, harvest_day; diurnal_config)
        @test reference.schedule_matches
        @test reference.yield > 0
        @info "Sub-daily weather attribution fixture" cft_id harvest_day production_yield = reference.yield
        flush(stderr)

        result = enzyme_weather_harvest_gradient(forcing, state, cft, parameters,
            climate, days, harvest_day; block_days = 4, diurnal_config)
        @test result.primal ≈ result.production_yield rtol = 1e-3 atol = 1e-5
        @test result.reverse_primal ≈ result.primal rtol = 1e-10 atol = 1e-12
        @test all(isfinite, result.gradient)
        @test forcing == original
        @test all(iszero, result.gradient[1:(first(days) - 1), :, :])
        @test all(iszero, result.gradient[harvest_day:end, :, :])
        @test any(!iszero, result.gradient[days, :, 1])

        # Temperature channel finite-difference cross-check: this is the
        # channel the sub-daily kernel actually reshapes (via `diurnal_range`
        # entering `diurnal_temperature` as a fixed Const auxiliary), so it is
        # the one most likely to expose a wrong `Enzyme.Const`/`Duplicated`
        # wiring if the sub-daily branch were differentiated incorrectly.
        direction = zeros(T, size(forcing))
        direction[first(days):(first(days) + 2), 1, 1] .= 1.0
        projection = sum(result.gradient .* direction)
        epsilon = T === Float32 ? T(0.05) : T(0.01)
        plus = weather_harvest_replay(forcing .+ epsilon .* direction, state,
            cft, parameters, climate, days, harvest_day; diurnal_config)
        minus = weather_harvest_replay(forcing .- epsilon .* direction, state,
            cft, parameters, climate, days, harvest_day; diurnal_config)
        @test plus.schedule_matches && minus.schedule_matches
        fd = (plus.yield - minus.yield) / (2epsilon)
        @info "Sub-daily weather directional validation" cft_id ad = projection fd
        flush(stderr)
        @test projection ≈ fd rtol = 0.03 atol = 2e-5

        # Guard rail: `diurnal_config` set but no `diurnal_range` on `climate`
        # must fail loudly and immediately, not silently ignore sub-daily.
        no_range_case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day)
        @test_throws ArgumentError weather_harvest_replay(no_range_case.forcing, no_range_case.state,
            no_range_case.cft, no_range_case.parameters, no_range_case.climate,
            no_range_case.days, no_range_case.harvest_day; diurnal_config)
    end
end
