using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop organ nitrogen redistribution" begin
    crop = init_crop(2, CuArray)
    crop.phenology.is_growing .= 1
    crop.nitrogen.total .= 0.7f0
    crop.carbon.leaf .= 2.0f0
    crop.carbon.root .= 3.0f0
    crop.carbon.storage .= 4.0f0
    crop.carbon.pool .= 1.0f0
    crop.nitrogen.leaf .= 10.0f0
    crop.nitrogen.root .= 20.0f0
    crop.nitrogen.storage .= 30.0f0
    crop.nitrogen.pool .= 40.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    first_sum = Array(crop.nitrogen.leaf .+ crop.nitrogen.root .+ crop.nitrogen.storage .+ crop.nitrogen.pool)
    first_leafn = Array(crop.nitrogen.leaf)

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    second_sum = Array(crop.nitrogen.leaf .+ crop.nitrogen.root .+ crop.nitrogen.storage .+ crop.nitrogen.pool)
    second_leafn = Array(crop.nitrogen.leaf)

    @test first_sum ≈ fill(0.7f0, 2) atol = 1.0f-6
    @test second_sum ≈ first_sum atol = 1.0f-6
    @test second_leafn ≈ first_leafn atol = 1.0f-7
end
