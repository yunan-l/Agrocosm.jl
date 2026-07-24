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
end
