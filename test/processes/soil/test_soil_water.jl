using Agrocosm
using Test

@testset "Staged daily soil-water update" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)

    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= Float32[40, 60, 100, 200, 200]
    soil.nitrogen.nitrate .= reshape(Float32[10, 20, 30, 40, 50], 5, 1)
    crop.water.interception .= 0.0f0
    crop.water.root_distribution .= 0.0f0
    crop.water.root_distribution[1, 1] = 1.0f0
    pedotransfer!(soil)

    storage_before = sum(soil.water.storage)
    nitrate_before = sum(soil.nitrogen.nitrate)
    water_availability_before = sum(soil.water.relative_content .* crop.water.root_distribution)
    soil_infiltration!(soil, crop, Float32[10.0])
    storage_after_infiltration = sum(soil.water.storage)
    water_availability_after = sum(soil.water.relative_content .* crop.water.root_distribution)

    @test storage_after_infiltration > storage_before
    @test storage_after_infiltration + soil.water.surface_runoff[1] +
          sum(soil.water.lateral_runoff) + soil.water.bottom_drainage[1] ≈ storage_before + 10.0f0 atol = 1.0f-4
    @test water_availability_after >= water_availability_before
    @test sum(soil.nitrogen.nitrate) + soil.nitrogen.leaching[1] ≈ nitrate_before atol = 1.0f-4

    crop.water.transpiration_layer .= 0.4f0
    soil.water.evaporation .= 0.2f0
    soil_evapotranspiration!(soil, crop)

    @test sum(soil.water.storage) ≈ storage_after_infiltration - 3.0f0 atol = 1.0f-5
end

@testset "Same-day rain pulse affects crop water stress" begin
    dry_soil = init_soil(1, soilparams.soildepth, identity)
    dry_crop = init_crop(1, identity)
    pet = init_pet(1, identity)

    dry_soil.properties.sand_fraction .= 0.4f0
    dry_soil.properties.clay_fraction .= 0.2f0
    dry_soil.water.storage .= Float32[25, 60, 100, 200, 200]
    dry_crop.phenology.is_growing .= 1
    dry_crop.carbon.root .= 50.0f0
    dry_crop.water.root_distribution .= 0.0f0
    dry_crop.water.root_distribution[1, 1] = 1.0f0
    dry_crop.water.canopy_wet .= 0.0f0
    dry_crop.water.interception .= 0.0f0
    pet.eeq .= 5.0f0
    pet.daylength .= 12.0f0
    pedotransfer!(dry_soil)

    wet_soil = deepcopy(dry_soil)
    wet_crop = deepcopy(dry_crop)
    soil_infiltration!(wet_soil, wet_crop, Float32[20.0])

    assimilation = Float32[5.0]
    co2 = Float32[40.0]
    transpiration!(assimilation, cft1, dry_crop, pet, dry_soil, co2)
    transpiration!(assimilation, cft1, wet_crop, pet, wet_soil, co2)

    @test wet_crop.water.supply_sum[1] > dry_crop.water.supply_sum[1]
    @test wet_crop.water.stress[1] > dry_crop.water.stress[1]
    @test sum(wet_crop.water.transpiration_layer) > sum(dry_crop.water.transpiration_layer)
end
