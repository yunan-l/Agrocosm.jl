using Agrocosm
using Test

function nitrogen_uptake_fixture(device = identity)
    crop = init_crop(1, device)
    soil = init_soil(1, soilparams.soildepth, device)

    crop.phenology.is_growing .= 1
    crop.nitrogen.total .= 0.1f0
    crop.carbon.leaf .= 20.0f0
    crop.carbon.root .= 100.0f0
    crop.nitrogen.leaf .= 0.0f0
    crop.nitrogen.root .= 0.0f0
    crop.nitrogen.demand_leaf .= 0.4f0
    crop.nitrogen.demand_total .= 1.0f0
    crop.water.root_distribution .= device(Float32[0.5, 0.3, 0.2, 0.0, 0.0])

    soil.water.relative_content .= 1.0f0
    soil.water.saturation_fraction .= 0.4f0
    soil.thermal.temperature .= 15.0f0
    soil.nitrogen.nitrate .= device(reshape(Float32[1.0, 0.0, 0.5, 0.0, 0.0], 5, 1))
    soil.nitrogen.ammonium .= device(reshape(Float32[0.0, 0.7, 0.0, 0.0, 0.0], 5, 1))

    return crop, soil
end

@testset "Separate NO3/NH4 uptake conserves nitrogen" begin
    crop, soil = nitrogen_uptake_fixture()
    plant_before = crop.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(crop, cft1, soil)

    plant_gain = crop.nitrogen.total[1] - plant_before
    soil_loss = soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium)

    @test crop.nitrogen.uptake[1] > 0.0f0
    @test crop.nitrogen.auto_fertilizer[1] == 0.0f0
    @test plant_gain ≈ crop.nitrogen.uptake[1] atol = 1.0f-6
    @test soil_loss ≈ crop.nitrogen.uptake[1] atol = 1.0f-6
    @test all(soil.nitrogen.nitrate .>= 0.0f0)
    @test all(soil.nitrogen.ammonium .>= 0.0f0)
    @test soil.nitrogen.nitrate[4, 1] == 0.0f0
    @test soil.nitrogen.ammonium[4, 1] == 0.0f0
end

@testset "Automatic fertilizer is an explicit external N input" begin
    crop, soil = nitrogen_uptake_fixture()
    plant_before = crop.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(crop, cft1, soil; auto_fertilizer = true)

    plant_gain = crop.nitrogen.total[1] - plant_before
    soil_loss = soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium)
    @test crop.nitrogen.total[1] ≈ crop.nitrogen.demand_total[1] atol = 1.0f-6
    @test crop.nitrogen.uptake[1] ≈ plant_gain atol = 1.0f-6
    @test plant_gain ≈ soil_loss + crop.nitrogen.auto_fertilizer[1] atol = 1.0f-6
    @test crop.nitrogen.auto_fertilizer[1] >= 0.0f0
    @test crop.nitrogen.stress[1] == 1.0f0
end

@testset "Nitrogen uptake respects remaining plant demand" begin
    crop, soil = nitrogen_uptake_fixture()
    crop.nitrogen.demand_total .= 0.105f0
    plant_before = crop.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(crop, cft1, soil)

    @test crop.nitrogen.uptake[1] ≈ 0.005f0 atol = 1.0f-6
    @test crop.nitrogen.total[1] - plant_before ≈ 0.005f0 atol = 1.0f-6
    @test soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium) ≈ 0.005f0 atol = 1.0f-6
end

@testset "Daily nitrogen uptake flux resets when demand is satisfied" begin
    crop, soil = nitrogen_uptake_fixture()
    crop.nitrogen.uptake .= 99.0f0
    crop.nitrogen.demand_total .= crop.nitrogen.total
    plant_before = crop.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(crop, cft1, soil)

    @test crop.nitrogen.uptake[1] == 0.0f0
    @test crop.nitrogen.total[1] == plant_before
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) == soil_before
end
