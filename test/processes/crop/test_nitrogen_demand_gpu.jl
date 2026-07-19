using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop nitrogen demand" begin
    crop, _, _, _ = init_crop(2, CuArray)
    crop.isgrowing .= 1
    crop.leafc .= 2.0f0
    crop.rootc .= 3.0f0
    crop.poolc .= 1.0f0
    crop.stoc .= 4.0f0

    ndemand_crop!(
        crop,
        cft1,
        CuArray(Float32[10.0, 5.0]),
        CuArray(Float32[25.0, 15.0]),
    )

    leaf_demand = Array(crop.ndemand_leaf)
    total_demand = Array(crop.ndemand_tot)
    @test all(isfinite, leaf_demand)
    @test all(total_demand .>= leaf_demand)
end
