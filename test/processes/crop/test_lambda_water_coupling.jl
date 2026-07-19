using Agrocosm
using Test

@testset "LPJ CO2 units and minimum conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop, _, _, _ = init_crop(1, identity)
    pet = init_pet(1, identity)

    soil.sand .= 0.4f0
    soil.clay .= 0.2f0
    soil.swc .= Float32[100, 150, 250, 400, 400]
    pedotransfer!(soil)

    crop.isgrowing .= 1
    crop.rootc .= 1000.0f0
    crop.rootdist .= 0.0f0
    crop.rootdist[1] = 1.0f0
    crop.fpar .= 0.8f0
    crop.canopy_wet .= 0.0f0
    pet.eeq .= 1.0f0
    pet.daylength .= 12.0f0

    adtmm = Float32[5.0]
    co2_pa = Float32[40.0]
    expected_gp = 1.6f0 * adtmm[1] /
                  (co2_pa[1] * 1.0f-5 * (1.0f0 - 0.8f0) * 12.0f0 * 3600.0f0) +
                  cft1.gmin * crop.fpar[1]

    transpiration!(adtmm, cft1, crop, pet, soil, co2_pa)

    @test crop.gp[1] ≈ expected_gp rtol = 1.0f-6
    @test crop.gp[1] > cft1.gmin * crop.fpar[1]
end

@testset "Inactive crop has zero conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop, _, _, _ = init_crop(1, identity)
    pet = init_pet(1, identity)

    crop.isgrowing .= 0
    crop.fpar .= 0.8f0
    pet.daylength .= 12.0f0
    pet.eeq .= 1.0f0

    transpiration!(Float32[5.0], cft1, crop, pet, soil, Float32[40.0])

    @test crop.gp[1] == 0.0f0
end
