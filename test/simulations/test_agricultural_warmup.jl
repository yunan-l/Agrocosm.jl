using Agrocosm
using Test

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
end
