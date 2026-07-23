using Agrocosm
using Test

@testset "Soil decomposition response" begin
    r = CropSoilDecompositionResponse(Float64)

    @testset "Lloyd-Taylor temperature response" begin
        ft(T) = Agrocosm.soil_decomposition_temperature_response(r, T)
        @test ft(10.0) ≈ 1.0                       # normalized to 1 at the 10 °C reference
        @test ft(20.0) > ft(10.0) > ft(0.0)        # increases with temperature
        @test ft(0.0) > 0                           # positive everywhere
        # closed form at 20 °C
        @test ft(20.0) ≈ exp(308.56 * (1 / 56.02 - 1 / (20.0 + 56.02 - 10.0))) rtol = 1e-12
    end

    @testset "moisture polynomial" begin
        fm(m) = Agrocosm.soil_decomposition_moisture_response(r, m)
        @test fm(0.0) ≈ 0.04021601
        @test fm(0.5) ≈ 0.04021601 + 0.71890122 * 0.5 + 4.26937932 * 0.25 + (-5.00505434) * 0.125 rtol = 1e-12
    end

    @testset "combined response bounded [0,1]" begin
        for T in (-10.0, 0.0, 10.0, 25.0, 40.0), m in (0.0, 0.25, 0.5, 0.75, 1.0)
            resp = Agrocosm.soil_decomposition_response(r, T, m)
            @test 0.0 ≤ resp ≤ 1.0
        end
    end
end
