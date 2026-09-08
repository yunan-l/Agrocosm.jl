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

        # G1/G2 in a real Enzyme pipeline: steps = 1 must degenerate exactly
        # onto the daily kernel, both in the ordinary replay and in the
        # gradient it feeds into `enzyme_weather_harvest_gradient`. This is
        # the same conservation property `test_diurnal.jl` checks in isolation,
        # now checked end-to-end through the actual AD path that had no
        # `diurnal_config` wiring at all before 2026-09-08.
        degenerate_config = DiurnalConfig(; steps = 1, shape = DIURNAL_SINUSOID)
        daily_case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day)
        degenerate_case = weather_attribution_fixture(cft_id; T, window_days = 8, sowing_day,
            diurnal_config = degenerate_config, diurnal_amplitude = T(10))
        daily_reference = weather_harvest_replay(daily_case.forcing, daily_case.state,
            daily_case.cft, daily_case.parameters, daily_case.climate, daily_case.days,
            daily_case.harvest_day)
        degenerate_reference = weather_harvest_replay(degenerate_case.forcing, degenerate_case.state,
            degenerate_case.cft, degenerate_case.parameters, degenerate_case.climate,
            degenerate_case.days, degenerate_case.harvest_day; diurnal_config = degenerate_config)
        @test degenerate_reference.yield == daily_reference.yield
        daily_gradient = enzyme_weather_harvest_gradient(daily_case.forcing, daily_case.state,
            daily_case.cft, daily_case.parameters, daily_case.climate, daily_case.days,
            daily_case.harvest_day; block_days = 4)
        degenerate_gradient = enzyme_weather_harvest_gradient(degenerate_case.forcing, degenerate_case.state,
            degenerate_case.cft, degenerate_case.parameters, degenerate_case.climate,
            degenerate_case.days, degenerate_case.harvest_day; block_days = 4, diurnal_config = degenerate_config)
        @test degenerate_gradient.gradient == daily_gradient.gradient

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
