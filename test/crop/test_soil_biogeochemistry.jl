using Agrocosm
using Test

@testset "Crop soil carbon biogeochemistry" begin
    bgc = CropSoilBiogeochemistry(Float64)
    tend(litter, fast, slow, response) = Agrocosm.soil_carbon_tendencies(bgc, litter, fast, slow, response)

    @testset "pools decompose; respiration positive" begin
        d_litter, d_fast, d_slow, het = tend(1.0, 5.0, 20.0, 1.0)
        @test d_litter < 0                       # litter is consumed
        @test het > 0                            # heterotrophic respiration
    end

    @testset "carbon is conserved (pool loss = respiration)" begin
        d_litter, d_fast, d_slow, het = tend(1.0, 5.0, 20.0, 1.0)
        # d(litter+fast+slow)/dt + heterotrophic respiration = 0
        @test d_litter + d_fast + d_slow + het ≈ 0.0 atol = 1e-18
    end

    @testset "no decomposition when the response is zero" begin
        @test tend(1.0, 5.0, 20.0, 0.0) == (0.0, 0.0, 0.0, 0.0)
    end

    @testset "warmer/wetter soil decomposes faster" begin
        _, _, _, het_cold = tend(1.0, 5.0, 20.0, Agrocosm.soil_decomposition_response(bgc.response, 5.0, 0.5))
        _, _, _, het_warm = tend(1.0, 5.0, 20.0, Agrocosm.soil_decomposition_response(bgc.response, 25.0, 0.5))
        @test het_warm > het_cold
    end

    @testset "mineral nitrogen transforms" begin
        ntend(nh4, no3, mineralization, T, wfps, orgc) =
            Agrocosm.soil_nitrogen_tendencies(bgc, nh4, no3, mineralization, T, wfps, orgc)

        # Mineralization feeds ammonium; when it is cold nitrification is off, so NH4 gains ≈ it.
        # (Denitrification retains a small constant temperature response below 0 °C, so NO3 still
        # declines slightly.)
        d_nh4, d_no3 = ntend(0.05, 0.05, 1.0e-8, -30.0, 0.6, 25.0)
        @test d_nh4 ≈ 1.0e-8 rtol = 1e-4
        @test d_no3 ≤ 0.0

        # Nitrification moves NH4 → NO3: with no mineralization, ammonium declines and nitrate rises.
        d_nh4b, d_no3b = ntend(0.1, 0.0, 0.0, 25.0, 0.6, 25.0)
        @test d_nh4b < 0
        @test d_no3b > 0

        # Denitrification removes nitrate (extra loss with organic carbon present, warm+wet).
        _, d_no3_lowC = ntend(0.0, 0.1, 0.0, 25.0, 0.9, 0.1)
        _, d_no3_highC = ntend(0.0, 0.1, 0.0, 25.0, 0.9, 50.0)
        @test d_no3_highC < d_no3_lowC   # more organic carbon → more denitrification loss
    end
end
