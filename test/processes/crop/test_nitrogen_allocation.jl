using Agrocosm
using Test

@testset "Crop organ nitrogen is a redistributed stock" begin
    crop, _, _, _ = init_crop(1, identity)
    crop.isgrowing .= 1
    crop.nitrogen .= 0.7f0
    crop.leafc .= 2.0f0
    crop.rootc .= 3.0f0
    crop.stoc .= 4.0f0
    crop.poolc .= 1.0f0

    # Deliberately inconsistent old organ values must not be accumulated.
    crop.leafn .= 10.0f0
    crop.rootn .= 20.0f0
    crop.ston .= 30.0f0
    crop.pooln .= 40.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    first_allocation = Float32[
        crop.leafn[1], crop.rootn[1], crop.ston[1], crop.pooln[1],
    ]

    @test sum(first_allocation) ≈ crop.nitrogen[1] atol = 1.0f-6
    @test all(first_allocation .>= 0.0f0)

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    second_allocation = Float32[
        crop.leafn[1], crop.rootn[1], crop.ston[1], crop.pooln[1],
    ]

    @test second_allocation ≈ first_allocation atol = 1.0f-7
    @test sum(second_allocation) ≈ crop.nitrogen[1] atol = 1.0f-6
end

@testset "Zero leaf carbon clears organ nitrogen safely" begin
    crop, _, _, _ = init_crop(1, identity)
    crop.isgrowing .= 1
    crop.nitrogen .= 0.7f0
    crop.leafc .= 0.0f0
    crop.rootc .= 3.0f0
    crop.leafn .= 1.0f0
    crop.rootn .= 1.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)

    @test crop.leafn[1] == 0.0f0
    @test crop.rootn[1] == 0.0f0
    @test crop.ston[1] == 0.0f0
    @test crop.pooln[1] == 0.0f0
end

@testset "Integrated automatic-fertilizer nitrogen cycle" begin
    crop, _, _, _ = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    crop.isgrowing .= 1
    crop.nitrogen .= 0.1f0
    crop.leafc .= 2.0f0
    crop.rootc .= 3.0f0
    crop.poolc .= 1.0f0
    crop.stoc .= 4.0f0

    crop_nitrogen!(
        crop,
        cft1,
        soil,
        Float32[10.0],
        Float32[25.0];
        auto_fertilizer = true,
    )

    organ_n = crop.leafn[1] + crop.rootn[1] + crop.pooln[1] + crop.ston[1]
    @test crop.nitrogen[1] ≈ crop.ndemand_tot[1] atol = 1.0f-6
    @test organ_n ≈ crop.nitrogen[1] atol = 1.0f-6
    @test crop.nautofertilizer[1] > 0.0f0
    @test crop.vscal[1] == 1.0f0
end
