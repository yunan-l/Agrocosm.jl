using Agrocosm
using Test

function warmup_array_snapshot(value)
    value isa AbstractArray && return Any[Array(value)]
    arrays = Any[]
    for name in fieldnames(typeof(value))
        append!(arrays, warmup_array_snapshot(getfield(value, name)))
    end
    return arrays
end

@testset "Agricultural warm-up" begin
    initial, _ = simulation_api_fixture(Float32)
    climate = climate_block(Float32, 365, 15, 1)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    initial_carbon = sum(Array(Agrocosm.soil_carbon_prognostic(simulation.state).litter))

    report = agricultural_warmup!(simulation, climate)

    @test report.years == 10
    @test report.days == 3650
    @test report.forcing_years == 1
    @test size(report.soil.total_carbon) == (10, 1)
    @test size(report.soil.total_nitrogen) == (10, 1)
    @test size(report.soil.fast_carbon) == (10, 1)
    @test size(report.soil.slow_nitrogen) == (10, 1)
    @test all(isfinite, report.soil.total_carbon)
    @test all(report.soil.total_carbon .>= 0)
    @test all(report.soil.total_nitrogen .>= 0)
    @test all(report.soil.mineral_nitrogen .>= 0)
    @test report.soil.litter_carbon[end] != initial_carbon
    @test simulation.simulated_days == 0
    @test size(simulation.output.crop.npp, 1) == 0
    @test all(iszero, Array(simulation.carbon_balance.residual))

    run_simulation!(simulation, climate; end_day = 3, spinup = false)
    @test simulation.simulated_days == 3
    @test size(simulation.output.crop.npp) == (3, 1)

    started = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    transition_day!(started, climate)
    @test_throws ArgumentError agricultural_warmup!(started, climate)
    @test_throws ArgumentError agricultural_warmup!(
        initialize_simulation(
            cft1, initial;
            indices = [1], T = Float32, days = 3,
            diagnostics = false, fertilizer = :no,
        ),
        climate_block(Float32, 364, 15, 1),
    )

    eager_climate = (
        temp = reshape(Float32[15 + 2 * (day > 365) for day in 1:730], 730, 1),
        prec = reshape(Float32[1 + (day % 11 == 0) for day in 1:730], 730, 1),
        sw = fill(180.0f0, 730, 1),
        lw = fill(-40.0f0, 730, 1),
        co2 = reshape(Float32[400 + (day > 365) for day in 1:730], 730, 1),
        co2_daily = true,
        backend_neutral = true,
    )
    ranges = [1:73, 74:289, 290:410, 411:619, 620:730]
    streamed_climate = [(
        temp = eager_climate.temp[range, :],
        prec = eager_climate.prec[range, :],
        sw = eager_climate.sw[range, :],
        lw = eager_climate.lw[range, :],
        co2 = eager_climate.co2[range, :],
        co2_daily = true,
        backend_neutral = true,
    ) for range in ranges]
    eager_simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    streamed_simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    eager_report = agricultural_warmup!(eager_simulation, eager_climate; years = 3)
    streamed_report = agricultural_warmup!(streamed_simulation, streamed_climate; years = 3)

    @test streamed_report == eager_report
    @test warmup_array_snapshot(streamed_simulation.state.prognostic) ==
        warmup_array_snapshot(eager_simulation.state.prognostic)
    @test streamed_simulation.simulated_days == 0
    @test size(streamed_simulation.output.crop.npp, 1) == 0
    @test all(iszero, Array(streamed_simulation.carbon_balance.residual))
    @test_throws ArgumentError agricultural_warmup!(
        initialize_simulation(
            cft1, initial;
            indices = [1], T = Float32, days = 3,
            diagnostics = false, fertilizer = :no,
        ),
        streamed_climate[1:2],
    )
end
