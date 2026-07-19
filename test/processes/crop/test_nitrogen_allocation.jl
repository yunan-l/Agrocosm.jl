using Agrocosm
using Test

@testset "Crop organ nitrogen is a redistributed stock" begin
    crop = init_crop(1, identity)
    crop.phenology.is_growing .= 1
    crop.nitrogen.total .= 0.7f0
    crop.carbon.leaf .= 2.0f0
    crop.carbon.root .= 3.0f0
    crop.carbon.storage .= 4.0f0
    crop.carbon.pool .= 1.0f0

    # Deliberately inconsistent old organ values must not be accumulated.
    crop.nitrogen.leaf .= 10.0f0
    crop.nitrogen.root .= 20.0f0
    crop.nitrogen.storage .= 30.0f0
    crop.nitrogen.pool .= 40.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    first_allocation = Float32[
        crop.nitrogen.leaf[1], crop.nitrogen.root[1], crop.nitrogen.storage[1], crop.nitrogen.pool[1],
    ]

    @test sum(first_allocation) ≈ crop.nitrogen.total[1] atol = 1.0f-6
    @test all(first_allocation .>= 0.0f0)

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)
    second_allocation = Float32[
        crop.nitrogen.leaf[1], crop.nitrogen.root[1], crop.nitrogen.storage[1], crop.nitrogen.pool[1],
    ]

    @test second_allocation ≈ first_allocation atol = 1.0f-7
    @test sum(second_allocation) ≈ crop.nitrogen.total[1] atol = 1.0f-6
end

@testset "Zero leaf carbon clears organ nitrogen safely" begin
    crop = init_crop(1, identity)
    crop.phenology.is_growing .= 1
    crop.nitrogen.total .= 0.7f0
    crop.carbon.leaf .= 0.0f0
    crop.carbon.root .= 3.0f0
    crop.nitrogen.leaf .= 1.0f0
    crop.nitrogen.root .= 1.0f0

    Agrocosm.allocate_crop_nitrogen!(crop, cft1)

    @test crop.nitrogen.leaf[1] == 0.0f0
    @test crop.nitrogen.root[1] == 0.0f0
    @test crop.nitrogen.storage[1] == 0.0f0
    @test crop.nitrogen.pool[1] == 0.0f0
end

@testset "Integrated automatic-fertilizer nitrogen cycle" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    crop.phenology.is_growing .= 1
    crop.nitrogen.total .= 0.1f0
    crop.carbon.leaf .= 2.0f0
    crop.carbon.root .= 3.0f0
    crop.carbon.pool .= 1.0f0
    crop.carbon.storage .= 4.0f0

    crop_nitrogen!(
        crop,
        cft1,
        soil,
        Float32[10.0],
        Float32[25.0];
        auto_fertilizer = true,
    )

    organ_n = crop.nitrogen.leaf[1] + crop.nitrogen.root[1] + crop.nitrogen.pool[1] + crop.nitrogen.storage[1]
    @test crop.nitrogen.total[1] ≈ crop.nitrogen.demand_total[1] atol = 1.0f-6
    @test organ_n ≈ crop.nitrogen.total[1] atol = 1.0f-6
    @test crop.nitrogen.auto_fertilizer[1] > 0.0f0
    @test crop.nitrogen.stress[1] == 1.0f0
end
