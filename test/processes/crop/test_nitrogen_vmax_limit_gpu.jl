using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop nitrogen limitation of Vmax" begin
    crop = init_crop(2, CuArray)
    crop.phenology.is_growing .= true
    crop.carbon.leaf .= 100.0f0
    crop.photosynthesis.potential_vmax .= CuArray(Float32[10, 10])
    crop.photosynthesis.vmax .= crop.photosynthesis.potential_vmax
    temperature = CuArray(Float32[25, 25])

    ndemand_crop!(crop, cft1, crop.photosynthesis.potential_vmax, temperature)
    @views crop.nitrogen.demand_leaf[2:2] .= cft1.ncleaf.low .* crop.carbon.leaf[2:2]
    limit_vmax_by_nitrogen!(crop, cft1, temperature)

    vmax = Array(crop.photosynthesis.vmax)
    limitation = Array(crop.photosynthesis.nitrogen_limitation)
    @test vmax[1] ≈ 10.0f0 atol = 2.0f-5
    @test 0.0f0 < vmax[2] < 10.0f0
    @test limitation[1] ≈ 1.0f0 atol = 2.0f-6
    @test limitation[2] < 1.0f-5
end
