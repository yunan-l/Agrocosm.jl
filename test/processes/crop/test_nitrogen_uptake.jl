using Agrocosm
using Test

function nitrogen_uptake_fixture(device = identity)
    crop, _, _, _ = init_crop(1, device)
    soil = init_soil(1, soilparams.soildepth, device)

    crop.isgrowing .= 1
    crop.nitrogen .= 0.1f0
    crop.leafc .= 20.0f0
    crop.rootc .= 100.0f0
    crop.leafn .= 0.0f0
    crop.rootn .= 0.0f0
    crop.ndemand_leaf .= 0.4f0
    crop.ndemand_tot .= 1.0f0
    crop.rootdist .= device(Float32[0.5, 0.3, 0.2, 0.0, 0.0])

    soil.w .= 1.0f0
    soil.wsat .= 0.4f0
    soil.temp .= 15.0f0
    soil.NO3 .= device(reshape(Float32[1.0, 0.0, 0.5, 0.0, 0.0], 5, 1))
    soil.NH4 .= device(reshape(Float32[0.0, 0.7, 0.0, 0.0, 0.0], 5, 1))

    return crop, soil
end

@testset "Separate NO3/NH4 uptake conserves nitrogen" begin
    crop, soil = nitrogen_uptake_fixture()
    plant_before = crop.nitrogen[1]
    soil_before = sum(soil.NO3) + sum(soil.NH4)

    nuptake_crop!(crop, cft1, soil)

    plant_gain = crop.nitrogen[1] - plant_before
    soil_loss = soil_before - sum(soil.NO3) - sum(soil.NH4)

    @test crop.nuptake[1] > 0.0f0
    @test plant_gain ≈ crop.nuptake[1] atol = 1.0f-6
    @test soil_loss ≈ crop.nuptake[1] atol = 1.0f-6
    @test all(soil.NO3 .>= 0.0f0)
    @test all(soil.NH4 .>= 0.0f0)
    @test soil.NO3[4, 1] == 0.0f0
    @test soil.NH4[4, 1] == 0.0f0
end

@testset "Nitrogen uptake respects remaining plant demand" begin
    crop, soil = nitrogen_uptake_fixture()
    crop.ndemand_tot .= 0.105f0
    plant_before = crop.nitrogen[1]
    soil_before = sum(soil.NO3) + sum(soil.NH4)

    nuptake_crop!(crop, cft1, soil)

    @test crop.nuptake[1] ≈ 0.005f0 atol = 1.0f-6
    @test crop.nitrogen[1] - plant_before ≈ 0.005f0 atol = 1.0f-6
    @test soil_before - sum(soil.NO3) - sum(soil.NH4) ≈ 0.005f0 atol = 1.0f-6
end

@testset "Daily nitrogen uptake flux resets when demand is satisfied" begin
    crop, soil = nitrogen_uptake_fixture()
    crop.nuptake .= 99.0f0
    crop.ndemand_tot .= crop.nitrogen
    plant_before = crop.nitrogen[1]
    soil_before = sum(soil.NO3) + sum(soil.NH4)

    nuptake_crop!(crop, cft1, soil)

    @test crop.nuptake[1] == 0.0f0
    @test crop.nitrogen[1] == plant_before
    @test sum(soil.NO3) + sum(soil.NH4) == soil_before
end
