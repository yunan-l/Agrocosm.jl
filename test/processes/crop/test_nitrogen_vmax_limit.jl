using Agrocosm
using Test

@testset "LPJmL crop nitrogen limitation of Vmax" begin
    crop = init_crop(2, identity)
    crop.state.phenology.is_growing .= true
    crop.state.carbon.leaf .= 100.0f0
    crop.auxiliary.photosynthesis.potential_vmax .= Float32[10, 10]
    crop.auxiliary.photosynthesis.vmax .= crop.auxiliary.photosynthesis.potential_vmax
    temperature = Float32[25, 25]

    # The demand equation and Vmax-limit equation are analytical inverses
    # when optimal leaf N is available.
    ndemand_crop!(crop, cft1, crop.auxiliary.photosynthesis.potential_vmax, temperature)
    optimal_leaf_n = copy(crop.auxiliary.stress.nitrogen_demand_leaf)
    crop.auxiliary.stress.nitrogen_demand_leaf[2] = cft1.ncleaf.low * crop.state.carbon.leaf[2]
    limit_vmax_by_nitrogen!(crop, cft1, temperature)

    @test crop.auxiliary.photosynthesis.vmax[1] ≈ 10.0f0 atol = 2.0f-5
    @test crop.auxiliary.photosynthesis.nitrogen_limitation[1] ≈ 1.0f0 atol = 2.0f-6
    @test 0.0f0 < crop.auxiliary.photosynthesis.vmax[2] < crop.auxiliary.photosynthesis.potential_vmax[2]
    @test crop.auxiliary.photosynthesis.nitrogen_limitation[2] < 1.0f-5
    @test crop.auxiliary.stress.nitrogen_demand_leaf[1] == optimal_leaf_n[1]
end
