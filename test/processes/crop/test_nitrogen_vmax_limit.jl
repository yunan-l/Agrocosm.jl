using Agrocosm
using Test

@testset "LPJmL crop nitrogen limitation of Vmax" begin
    crop = init_crop(2, identity)
    crop.phenology.is_growing .= true
    crop.carbon.leaf .= 100.0f0
    crop.photosynthesis.potential_vmax .= Float32[10, 10]
    crop.photosynthesis.vmax .= crop.photosynthesis.potential_vmax
    temperature = Float32[25, 25]

    # The demand equation and Vmax-limit equation are analytical inverses
    # when optimal leaf N is available.
    ndemand_crop!(crop, cft1, crop.photosynthesis.potential_vmax, temperature)
    optimal_leaf_n = copy(crop.nitrogen.demand_leaf)
    crop.nitrogen.demand_leaf[2] = cft1.ncleaf.low * crop.carbon.leaf[2]
    limit_vmax_by_nitrogen!(crop, cft1, temperature)

    @test crop.photosynthesis.vmax[1] ≈ 10.0f0 atol = 2.0f-5
    @test crop.photosynthesis.nitrogen_limitation[1] ≈ 1.0f0 atol = 2.0f-6
    @test 0.0f0 < crop.photosynthesis.vmax[2] < crop.photosynthesis.potential_vmax[2]
    @test crop.photosynthesis.nitrogen_limitation[2] < 1.0f-5
    @test crop.nitrogen.demand_leaf[1] == optimal_leaf_n[1]
end
