using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop organ nitrogen redistribution" begin
    crop = init_crop(2, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.7f0
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.storage .= 4.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.nitrogen.leaf .= 10.0f0
    crop.state.nitrogen.root .= 20.0f0
    crop.state.nitrogen.storage .= 30.0f0
    crop.state.nitrogen.pool .= 40.0f0

    Agrocosm.allocate_crop_nitrogen!(state, cft1)
    first_sum = Array(crop.state.nitrogen.leaf .+ crop.state.nitrogen.root .+ crop.state.nitrogen.storage .+ crop.state.nitrogen.pool)
    first_leafn = Array(crop.state.nitrogen.leaf)

    Agrocosm.allocate_crop_nitrogen!(state, cft1)
    second_sum = Array(crop.state.nitrogen.leaf .+ crop.state.nitrogen.root .+ crop.state.nitrogen.storage .+ crop.state.nitrogen.pool)
    second_leafn = Array(crop.state.nitrogen.leaf)

    @test first_sum ≈ fill(0.7f0, 2) atol = 1.0f-6
    @test second_sum ≈ first_sum atol = 1.0f-6
    @test second_leafn ≈ first_leafn atol = 1.0f-7
end

@testset "CUDA living leafless crops retain organ nitrogen" begin
    crop = init_crop(4, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= CuArray(Int32[1, 1, 1, 0])
    crop.state.nitrogen.total .= 2.0f0
    crop.state.carbon.leaf .= CuArray(Float32[0, 1e-8, 1e-7, 0])
    crop.state.carbon.root .= 100.0f0
    crop.state.carbon.storage .= 200.0f0
    crop.state.carbon.pool .= 100.0f0
    crop.state.nitrogen.leaf .= 0.0f0
    crop.state.nitrogen.root .= 0.5f0
    crop.state.nitrogen.storage .= 1.0f0
    crop.state.nitrogen.pool .= 0.5f0

    Agrocosm.allocate_crop_nitrogen!(state, cft1)

    @test all(iszero, Array(crop.state.nitrogen.leaf))
    @test Array(crop.state.nitrogen.root) == Float32[0.5, 0.5, 0.5, 0]
    @test Array(crop.state.nitrogen.storage) == Float32[1, 1, 1, 0]
    @test Array(crop.state.nitrogen.pool) == Float32[0.5, 0.5, 0.5, 0]
end
