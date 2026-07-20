using Agrocosm
using Test

@testset "Harvest carbon and nitrogen routing is conservative" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    residue_fraction = Float32[0.6]
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]

    crop.carbon.leaf .= 2.0f0
    crop.carbon.root .= 4.0f0
    crop.carbon.storage .= 5.0f0
    crop.carbon.pool .= 3.0f0
    crop.carbon.biomass .= 14.0f0
    crop.nitrogen.leaf .= 0.2f0
    crop.nitrogen.root .= 0.4f0
    crop.nitrogen.storage .= 0.5f0
    crop.nitrogen.pool .= 0.3f0
    crop.nitrogen.total .= 1.4f0
    crop.phenology.is_growing .= Int32(1)
    crop.phenology.harvesting_previous .= false
    crop.phenology.harvesting .= true

    carbon_before = crop.carbon.leaf[1] + crop.carbon.root[1] +
        crop.carbon.storage[1] + crop.carbon.pool[1]
    nitrogen_before = crop.nitrogen.total[1]
    harvest_crop!(crop.calendar, crop, soil, output, residue_fraction, 150)

    carbon_export = crop.carbon.yield[1] +
        (crop.carbon.leaf[1] + crop.carbon.pool[1]) * (1 - residue_fraction[1])
    nitrogen_export = crop.nitrogen.harvest_export[1]
    carbon_residue = sum(soil.carbon.input)
    nitrogen_residue = sum(soil.nitrogen.input)
    @test crop.calendar.harvest_callback[1] == 1
    @test crop.carbon.yield[1] == 5.0f0
    @test carbon_export == 7.0f0
    @test carbon_residue == 7.0f0
    @test carbon_export + carbon_residue == carbon_before
    @test nitrogen_export == 0.7f0
    @test nitrogen_residue ≈ 0.7f0 atol = 1.0f-7
    @test nitrogen_export + nitrogen_residue ≈ nitrogen_before atol = 2.0f-7

    Agrocosm.route_harvest_residues!(soil, crop.calendar)
    @test sum(soil.carbon.litter) == carbon_residue
    @test sum(soil.nitrogen.litter) == nitrogen_residue

    # The production loop clears the harvested persistent crop object later on
    # the same day through the inactive carbon and nitrogen kernel branches.
    carbon_allocation!(cft1, crop, crop.photosynthesis)
    crop_nitrogen!(
        crop, cft1, soil, zeros(Float32, 1), Float32[15];
        auto_fertilizer = false,
    )
    @test crop.carbon.biomass[1] == 0.0f0
    @test crop.carbon.leaf[1] == 0.0f0
    @test crop.carbon.root[1] == 0.0f0
    @test crop.carbon.storage[1] == 0.0f0
    @test crop.carbon.pool[1] == 0.0f0
    @test crop.nitrogen.total[1] == 0.0f0
end
