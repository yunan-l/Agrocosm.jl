using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA separate NO3/NH4 uptake conserves nitrogen" begin
    crop, _, _, _ = init_crop(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)

    crop.isgrowing .= 1
    crop.nitrogen .= 0.1f0
    crop.leafc .= 20.0f0
    crop.rootc .= 100.0f0
    crop.ndemand_leaf .= 0.4f0
    crop.ndemand_tot .= 1.0f0
    crop.rootdist .= CuArray(Float32[0.5, 0.3, 0.2, 0.0, 0.0])
    soil.w .= 1.0f0
    soil.wsat .= 0.4f0
    soil.temp .= 15.0f0
    soil.NO3 .= CuArray(reshape(Float32[1.0, 0.0, 0.5, 0.0, 0.0], 5, 1))
    soil.NH4 .= CuArray(reshape(Float32[0.0, 0.7, 0.0, 0.0, 0.0], 5, 1))

    plant_before = Array(crop.nitrogen)[1]
    soil_before = sum(Array(soil.NO3)) + sum(Array(soil.NH4))
    nuptake_crop!(crop, cft1, soil)

    uptake = Array(crop.nuptake)[1]
    plant_gain = Array(crop.nitrogen)[1] - plant_before
    soil_loss = soil_before - sum(Array(soil.NO3)) - sum(Array(soil.NH4))

    @test uptake > 0.0f0
    @test plant_gain ≈ uptake atol = 1.0f-6
    @test soil_loss ≈ uptake atol = 1.0f-6
    @test all(Array(soil.NO3) .>= 0.0f0)
    @test all(Array(soil.NH4) .>= 0.0f0)
end
