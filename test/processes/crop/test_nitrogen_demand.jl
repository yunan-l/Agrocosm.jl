using Agrocosm
using Test

@testset "Crop nitrogen demand follows LPJmL ndemand_crop" begin
    crop, _, _, _ = init_crop(1, identity)
    crop.isgrowing .= 1
    crop.leafc .= 2.0f0
    crop.rootc .= 3.0f0
    crop.poolc .= 1.0f0
    crop.stoc .= 4.0f0
    vmax = Float32[10.0]
    temp = Float32[25.0]

    ndemand_crop!(crop, cft1, vmax, temp)

    expected_leaf = lpjmlparams.p * 1.0f-3 * vmax[1] /
                    (86400.0f0 * 12.0f0 * 1.0f-6) +
                    cft1.ncleaf.low * crop.leafc[1]
    expected_nc = clamp(
        expected_leaf / crop.leafc[1],
        cft1.ncleaf.low,
        cft1.ncleaf.high,
    )
    expected_total = expected_leaf + expected_nc * (
        crop.rootc[1] / cft1.ratio.root +
        crop.poolc[1] / cft1.ratio.pool +
        crop.stoc[1] / cft1.ratio.sto
    )

    @test crop.ndemand_leaf[1] ≈ expected_leaf atol = 1.0f-6
    @test crop.ndemand_tot[1] ≈ expected_total atol = 1.0f-6
    @test crop.ndemand_tot[1] >= crop.ndemand_leaf[1]
end

@testset "Inactive crop has no nitrogen demand" begin
    crop, _, _, _ = init_crop(1, identity)
    crop.ndemand_leaf .= 9.0f0
    crop.ndemand_tot .= 9.0f0

    ndemand_crop!(crop, cft1, Float32[10.0], Float32[25.0])

    @test crop.ndemand_leaf[1] == 0.0f0
    @test crop.ndemand_tot[1] == 0.0f0
end
