using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop nitrogen demand" begin
    crop = init_crop(2, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.carbon.storage .= 4.0f0
    crop.auxiliary.photosynthesis.lambda .= CuArray(Float32[0.8, 0.0])

    ndemand_crop!(
        state,
        cft1,
        CuArray(Float32[10.0, 5.0]),
        CuArray(Float32[25.0, 15.0]),
        require_active_photosynthesis = true,
    )

    leaf_demand = Array(crop.auxiliary.stress.nitrogen_demand_leaf)
    total_demand = Array(crop.auxiliary.stress.nitrogen_demand_total)
    @test all(isfinite, leaf_demand)
    @test total_demand[1] >= leaf_demand[1] > 0.0f0
    @test total_demand[2] == leaf_demand[2] == 0.0f0
end
