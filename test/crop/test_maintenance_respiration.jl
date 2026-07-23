using Agrocosm
using Test

@testset "Crop maintenance respiration" begin
    m = CropMaintenanceRespiration(Float64)   # respcoeff=0.8, k=0.0548, e0=308.56, temp_response=56.02

    @testset "Lloyd-Taylor temperature response" begin
        @test Agrocosm.maintenance_temperature_response(m, 10.0) ≈ 1.0        # normalized to 1 at 10 °C
        @test Agrocosm.maintenance_temperature_response(m, -20.0) == 0.0      # gated below -15 °C
        # monotone increasing with temperature over the active range
        @test Agrocosm.maintenance_temperature_response(m, 25.0) > Agrocosm.maintenance_temperature_response(m, 5.0)
    end

    @testset "organ and total maintenance respiration" begin
        # at 10 °C the temperature response is 1, so respiration is carbon·respcoeff·k·nc_ratio
        g = 1.0
        @test Agrocosm.organ_maintenance_respiration(m, 50.0, 1 / 30, g) ≈ 50.0 * 0.8 * 0.0548 * (1 / 30)
        total = Agrocosm.maintenance_respiration(m, 50.0, 5.0, 10.0, 10.0, 10.0)
        expected =
            50.0 * 0.8 * 0.0548 * (1 / 30) +
            5.0 * 0.8 * 0.0548 * (1 / 100) +
            10.0 * 0.8 * 0.0548 * (1 / 100)
        @test total ≈ expected rtol = 1e-9
        @test total ≥ 0
    end

    @testset "warmer soil/air raise respiration; zero carbon gives zero" begin
        cold = Agrocosm.maintenance_respiration(m, 50.0, 5.0, 10.0, 5.0, 5.0)
        warm = Agrocosm.maintenance_respiration(m, 50.0, 5.0, 10.0, 25.0, 25.0)
        @test warm > cold
        @test Agrocosm.maintenance_respiration(m, 0.0, 0.0, 0.0, 20.0, 20.0) == 0.0
    end
end
