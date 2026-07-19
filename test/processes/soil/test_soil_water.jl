using Agrocosm
using Test

@testset "Staged daily soil-water update" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop, _, _, _ = init_crop(1, identity)

    soil.sand .= 0.4f0
    soil.clay .= 0.2f0
    soil.swc .= Float32[40, 60, 100, 200, 200]
    crop.intercep .= 0.0f0
    crop.rootdist .= 0.0f0
    crop.rootdist[1, 1] = 1.0f0
    pedotransfer!(soil)

    storage_before = sum(soil.swc)
    water_availability_before = sum(soil.w .* crop.rootdist)
    soil_infiltration!(soil, crop, Float32[10.0])
    storage_after_infiltration = sum(soil.swc)
    water_availability_after = sum(soil.w .* crop.rootdist)

    @test storage_after_infiltration > storage_before
    @test storage_after_infiltration + soil.srunoff[1] +
          sum(soil.lrunoff) + soil.outflux_f[1] ≈ storage_before + 10.0f0 atol = 1.0f-4
    @test water_availability_after >= water_availability_before

    crop.trans_layer .= 0.4f0
    soil.evap .= 0.2f0
    soil_evapotranspiration!(soil, crop)

    @test sum(soil.swc) ≈ storage_after_infiltration - 3.0f0 atol = 1.0f-5
end

@testset "Same-day rain pulse affects crop water stress" begin
    dry_soil = init_soil(1, soilparams.soildepth, identity)
    dry_crop, _, _, _ = init_crop(1, identity)
    pet = init_pet(1, identity)

    dry_soil.sand .= 0.4f0
    dry_soil.clay .= 0.2f0
    dry_soil.swc .= Float32[25, 60, 100, 200, 200]
    dry_crop.isgrowing .= 1
    dry_crop.rootc .= 50.0f0
    dry_crop.rootdist .= 0.0f0
    dry_crop.rootdist[1, 1] = 1.0f0
    dry_crop.canopy_wet .= 0.0f0
    dry_crop.intercep .= 0.0f0
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

    @test wet_crop.w_supplysum[1] > dry_crop.w_supplysum[1]
    @test wet_crop.wscal[1] > dry_crop.wscal[1]
    @test sum(wet_crop.trans_layer) > sum(dry_crop.trans_layer)
end
