using Agrocosm
using Test

@testset "Crop growth respiration and NPP" begin
    g = CropGrowthRespiration(Float64)   # r_growth = 0.25
    Rg(gpp, rm) = Agrocosm.growth_respiration(g, gpp, rm)
    NPP(gpp, rm) = Agrocosm.net_primary_production(g, gpp, rm)

    @testset "positive net assimilate" begin
        @test Rg(10.0, 2.0) ≈ 0.25 * 8.0           # r_growth·(gross − maintenance)
        @test NPP(10.0, 2.0) ≈ 0.75 * 8.0          # (1 − r_growth)·(gross − maintenance)
        @test NPP(10.0, 2.0) ≈ 10.0 - 2.0 - Rg(10.0, 2.0)   # NPP = GPP − Rm − Rg
    end

    @testset "no growth respiration when maintenance exceeds assimilation" begin
        @test Rg(1.0, 3.0) == 0.0                  # gross < maintenance → no growth respiration
        @test NPP(1.0, 3.0) ≈ 1.0 - 3.0            # net loss = gross − maintenance
    end

    @testset "parameter sensitivity" begin
        g2 = CropGrowthRespiration(Float64; r_growth = 0.4)
        @test Agrocosm.growth_respiration(g2, 10.0, 2.0) ≈ 0.4 * 8.0
        @test Agrocosm.net_primary_production(g2, 10.0, 2.0) ≈ 0.6 * 8.0
    end
end
