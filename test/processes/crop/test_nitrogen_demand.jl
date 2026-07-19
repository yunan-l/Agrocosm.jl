using Agrocosm
using Test

@testset "Crop nitrogen demand follows LPJmL ndemand_crop" begin
    crop = init_crop(1, identity)
    crop.phenology.is_growing .= 1
    crop.carbon.leaf .= 2.0f0
    crop.carbon.root .= 3.0f0
    crop.carbon.pool .= 1.0f0
    crop.carbon.storage .= 4.0f0
    vmax = Float32[10.0]
    temp = Float32[25.0]

    ndemand_crop!(crop, cft1, vmax, temp)

    expected_leaf = lpjmlparams.p * 1.0f-3 * vmax[1] /
                    (86400.0f0 * 12.0f0 * 1.0f-6) +
                    cft1.ncleaf.low * crop.carbon.leaf[1]
    expected_nc = clamp(
        expected_leaf / crop.carbon.leaf[1],
        cft1.ncleaf.low,
        cft1.ncleaf.high,
    )
    expected_total = expected_leaf + expected_nc * (
        crop.carbon.root[1] / cft1.ratio.root +
        crop.carbon.pool[1] / cft1.ratio.pool +
        crop.carbon.storage[1] / cft1.ratio.sto
    )

    @test crop.nitrogen.demand_leaf[1] ≈ expected_leaf atol = 1.0f-6
    @test crop.nitrogen.demand_total[1] ≈ expected_total atol = 1.0f-6
    @test crop.nitrogen.demand_total[1] >= crop.nitrogen.demand_leaf[1]
end

@testset "Inactive crop has no nitrogen demand" begin
    crop = init_crop(1, identity)
    crop.nitrogen.demand_leaf .= 9.0f0
    crop.nitrogen.demand_total .= 9.0f0

    ndemand_crop!(crop, cft1, Float32[10.0], Float32[25.0])

    @test crop.nitrogen.demand_leaf[1] == 0.0f0
    @test crop.nitrogen.demand_total[1] == 0.0f0
end
