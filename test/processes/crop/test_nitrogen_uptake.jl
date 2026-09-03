using Agrocosm
using Test

function nitrogen_uptake_fixture(device = identity)
    crop = init_crop(1, device)
    soil = init_soil(1, soilparams.soildepth, device)

    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.1f0
    crop.state.carbon.leaf .= 20.0f0
    crop.state.carbon.root .= 100.0f0
    crop.state.nitrogen.leaf .= 0.0f0
    crop.state.nitrogen.root .= 0.0f0
    crop.auxiliary.stress.nitrogen_demand_leaf .= 0.4f0
    crop.auxiliary.stress.nitrogen_demand_total .= 1.0f0
    crop.auxiliary.root.distribution .= device(Float32[0.5, 0.3, 0.2, 0.0, 0.0])

    soil.water.relative_content .= 1.0f0
    soil.water.saturation_fraction .= 0.4f0
    soil.thermal.temperature .= 15.0f0
    soil.nitrogen.nitrate .= device(reshape(Float32[1.0, 0.0, 0.5, 0.0, 0.0], 5, 1))
    soil.nitrogen.ammonium .= device(reshape(Float32[0.0, 0.7, 0.0, 0.0, 0.0], 5, 1))

    return crop, soil, test_model_state(crop, soil)
end

@testset "LPJmL photosynthesis gate skips daily nitrogen acquisition" begin
    crop, soil, state = nitrogen_uptake_fixture()
    crop.auxiliary.photosynthesis.lambda .= 0.0f0
    crop.state.nitrogen.sufficiency .= 0.37f0
    crop.fluxes.nitrogen.uptake .= 99.0f0
    plant_before = copy(crop.state.nitrogen.total)
    nitrate_before = copy(soil.nitrogen.nitrate)
    ammonium_before = copy(soil.nitrogen.ammonium)

    nuptake_crop!(state, cft1, state; require_active_photosynthesis = true)

    @test crop.state.nitrogen.total == plant_before
    @test soil.nitrogen.nitrate == nitrate_before
    @test soil.nitrogen.ammonium == ammonium_before
    @test crop.fluxes.nitrogen.uptake[1] == 0.0f0
    @test crop.state.nitrogen.sufficiency[1] == 0.37f0
end

@testset "Separate NO3/NH4 uptake conserves nitrogen" begin
    crop, soil, state = nitrogen_uptake_fixture()
    plant_before = crop.state.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(state, cft1, state)

    plant_gain = crop.state.nitrogen.total[1] - plant_before
    soil_loss = soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium)

    @test crop.fluxes.nitrogen.uptake[1] > 0.0f0
    @test crop.fluxes.nitrogen.auto_fertilizer[1] == 0.0f0
    @test plant_gain ≈ crop.fluxes.nitrogen.uptake[1] atol = 1.0f-6
    @test soil_loss ≈ crop.fluxes.nitrogen.uptake[1] atol = 1.0f-6
    @test all(soil.nitrogen.nitrate .>= 0.0f0)
    @test all(soil.nitrogen.ammonium .>= 0.0f0)
    @test soil.nitrogen.nitrate[4, 1] == 0.0f0
    @test soil.nitrogen.ammonium[4, 1] == 0.0f0
end

@testset "Automatic fertilizer is an explicit external N input" begin
    crop, soil, state = nitrogen_uptake_fixture()
    plant_before = crop.state.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(state, cft1, state; auto_fertilizer = true)

    plant_gain = crop.state.nitrogen.total[1] - plant_before
    soil_loss = soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium)
    @test crop.state.nitrogen.total[1] ≈ crop.auxiliary.stress.nitrogen_demand_total[1] atol = 1.0f-6
    @test crop.fluxes.nitrogen.uptake[1] ≈ plant_gain atol = 1.0f-6
    @test plant_gain ≈ soil_loss + crop.fluxes.nitrogen.auto_fertilizer[1] atol = 1.0f-6
    @test crop.fluxes.nitrogen.auto_fertilizer[1] >= 0.0f0
    @test crop.state.nitrogen.sufficiency[1] == 1.0f0
end

@testset "Nitrogen uptake respects remaining plant demand" begin
    crop, soil, state = nitrogen_uptake_fixture()
    crop.auxiliary.stress.nitrogen_demand_total .= 0.105f0
    plant_before = crop.state.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(state, cft1, state)

    @test crop.fluxes.nitrogen.uptake[1] ≈ 0.005f0 atol = 1.0f-6
    @test crop.state.nitrogen.total[1] - plant_before ≈ 0.005f0 atol = 1.0f-6
    @test soil_before - sum(soil.nitrogen.nitrate) - sum(soil.nitrogen.ammonium) ≈ 0.005f0 atol = 1.0f-6
end

@testset "Daily nitrogen uptake flux resets when demand is satisfied" begin
    crop, soil, state = nitrogen_uptake_fixture()
    crop.fluxes.nitrogen.uptake .= 99.0f0
    crop.auxiliary.stress.nitrogen_demand_total .= crop.state.nitrogen.total
    plant_before = crop.state.nitrogen.total[1]
    soil_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    nuptake_crop!(state, cft1, state)

    @test crop.fluxes.nitrogen.uptake[1] == 0.0f0
    @test crop.state.nitrogen.total[1] == plant_before
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) == soil_before
end

@testset "Soybean biological fixation follows LPJmL NPP and soil controls" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.1f0
    crop.state.nitrogen.leaf .= 0.04f0
    crop.state.nitrogen.root .= 0.03f0
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.auxiliary.stress.nitrogen_demand_leaf .= 0.4f0
    crop.auxiliary.stress.nitrogen_demand_total .= 0.6f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1:2] .= 0.5f0
    crop.fluxes.carbon.net_assimilation .= 12.0f0
    soil.water.relative_content .= 0.0f0
    soil.water.relative_content[1:2] .= 0.8f0
    soil.thermal.temperature .= 25.0f0

    nuptake_crop!(
        state, cft9, state;
        auto_fertilizer = true,
        biological_fixation = true,
    )

    @test crop.fluxes.nitrogen.biological_fixation[1] ≈ 0.5f0
    @test crop.fluxes.carbon.biological_fixation_cost[1] ≈ 3.0f0
    @test crop.state.nitrogen.total[1] ≈ 0.6f0
    @test crop.fluxes.nitrogen.uptake[1] ≈ 0.5f0
    @test crop.fluxes.nitrogen.auto_fertilizer[1] == 0.0f0
    @test crop.state.nitrogen.sufficiency[1] == 1.0f0
end

@testset "Biological fixation responses are bounded" begin
    @test Agrocosm.compute_bnf_temperature_response(
        5.0f0, 5.0f0, 20.0f0, 35.0f0, 44.0f0,
    ) == 0.0f0
    @test Agrocosm.compute_bnf_temperature_response(
        27.0f0, 5.0f0, 20.0f0, 35.0f0, 44.0f0,
    ) == 1.0f0
    @test Agrocosm.compute_bnf_water_response(0.2f0, 0.2f0, 0.8f0) == 0.0f0
    @test Agrocosm.compute_bnf_water_response(0.8f0, 0.2f0, 0.8f0) == 1.0f0
end
