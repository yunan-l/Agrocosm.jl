using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop organ nitrogen redistribution" begin
    crop, _, _, _ = init_crop(2, CuArray)
    crop.isgrowing .= 1
    crop.nitrogen .= 0.7f0
    crop.leafc .= 2.0f0
    crop.rootc .= 3.0f0
    crop.stoc .= 4.0f0
    crop.poolc .= 1.0f0
    crop.leafn .= 10.0f0
    crop.rootn .= 20.0f0
    crop.ston .= 30.0f0
    crop.pooln .= 40.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    first_sum = Array(crop.leafn .+ crop.rootn .+ crop.ston .+ crop.pooln)
    first_leafn = Array(crop.leafn)

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    second_sum = Array(crop.leafn .+ crop.rootn .+ crop.ston .+ crop.pooln)
    second_leafn = Array(crop.leafn)

    @test first_sum ≈ fill(0.7f0, 2) atol = 1.0f-6
    @test second_sum ≈ first_sum atol = 1.0f-6
    @test second_leafn ≈ first_leafn atol = 1.0f-7
end
