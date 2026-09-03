using Agrocosm
using Test

@testset "Surface-litter water capacity and interception" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(soil)
    litter_carbon = 20.0f0 * 0.42f0 * 71.1f0
    soil.carbon.litter[1, 1] = litter_carbon

    update_surface_litter_properties!(state)

    expected_dry_matter = litter_carbon / 0.42f0
    expected_capacity = 2.0f-3 * expected_dry_matter
    expected_cover = 1.0f0 - exp(-6.0f-3 * expected_dry_matter)
    @test soil.surface_litter.depth[1] ≈ 0.02f0 atol = 1.0f-6
    @test soil.surface_litter.water_capacity[1] ≈ expected_capacity
    @test soil.surface_litter.cover[1] ≈ expected_cover

    soil.water.infiltration .= 10.0f0
    surface_litter_interception!(state)
    @test soil.surface_litter.interception[1] ≈ expected_capacity
    @test soil.surface_litter.water_storage[1] ≈ expected_capacity
    @test soil.water.infiltration[1] ≈ 10.0f0 - expected_capacity
end

@testset "Reduced litter capacity conserves water" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.carbon.litter[1, 1] = 100.0f0
    update_surface_litter_properties!(state)
    soil.surface_litter.water_storage .= soil.surface_litter.water_capacity
    soil.water.storage[1, 1] = 5.0f0
    water_before = soil.water.storage[1, 1] +
                   soil.surface_litter.water_storage[1]

    soil.carbon.litter[1, 1] = 50.0f0
    update_surface_litter_properties!(state)

    water_after = soil.water.storage[1, 1] +
                  soil.surface_litter.water_storage[1]
    @test water_after ≈ water_before atol = 1.0f-6
    @test soil.surface_litter.water_storage[1] ≈
          soil.surface_litter.water_capacity[1]
end

@testset "Wet litter evaporation" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    soil.carbon.litter[1, 1] = 20.0f0 * 0.42f0 * 71.1f0
    update_surface_litter_properties!(state)
    soil.surface_litter.water_storage .= soil.surface_litter.water_capacity
    storage_before = soil.surface_litter.water_storage[1]

    soil.water.storage .= 50.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    evaporation!(Float32[2.0], state, state)

    @test soil.surface_litter.evaporation[1] > 0.0f0
    @test soil.surface_litter.water_storage[1] ≈
          storage_before - soil.surface_litter.evaporation[1]
    @test soil.surface_litter.water_storage[1] >= 0.0f0
end

@testset "LPJmL managed-land soil evaporation" begin
    base_soil = init_soil(3, soilparams.soildepth, identity)
    base_crop = init_crop(3, identity)

    base_soil.properties.layer_depth .= Float32[200, 300, 500, 1000, 1000]
    base_soil.water.relative_content .= 0.0f0
    base_soil.water.relative_content[1:2, :] .= Float32[
        0.25 0.25 1.0f-5
        0.25 0.25 1.0f-5
    ]
    base_soil.water.free_water .= 0.0f0
    base_soil.water.holding_capacity_storage .= 100.0f0
    base_soil.water.ice_storage .= 0.0f0
    base_soil.water.ice_storage[1, 3] = 50.0f0
    base_soil.water.available_ice_storage .= 0.0f0
    base_soil.water.available_ice_storage[1, 3] = 50.0f0
    base_soil.water.free_ice_storage .= 0.0f0
    base_soil.surface_litter.cover .= Float32[0.0, 1.0, 0.0]
    base_soil.surface_litter.water_capacity .= 0.0f0
    base_soil.surface_litter.water_storage .= 0.0f0

    base_crop.auxiliary.canopy.fpar .= Float32[1.0, 0.0, 0.0]
    base_crop.auxiliary.canopy.canopy_wet .= 0.0f0
    base_crop.fluxes.water.transpiration_layer .= 0.0f0
    base_crop.fluxes.water.transpiration_layer[1:2, 3] .= 2.5f-4

    default_soil, default_crop = deepcopy.((base_soil, base_crop))
    explicit_true_soil, explicit_true_crop = deepcopy.((base_soil, base_crop))
    legacy_soil, legacy_crop = deepcopy.((base_soil, base_crop))
    default_state = test_model_state(default_crop, default_soil)
    explicit_true_state = test_model_state(explicit_true_crop, explicit_true_soil)
    legacy_state = test_model_state(legacy_crop, legacy_soil)
    equilibrium_evaporation = fill(2.0f0, 3)

    evaporation!(equilibrium_evaporation, default_state, default_state)
    evaporation!(
        equilibrium_evaporation, explicit_true_state, explicit_true_state;
        lpjml_managed_evaporation = true,
    )
    evaporation!(
        equilibrium_evaporation, legacy_state, legacy_state;
        lpjml_managed_evaporation = false,
    )

    @test reinterpret(UInt32, vec(default_soil.water.evaporation)) ==
          reinterpret(UInt32, vec(explicit_true_soil.water.evaporation))
    @test reinterpret(UInt32, default_soil.surface_litter.evaporation) ==
          reinterpret(UInt32, explicit_true_soil.surface_litter.evaporation)

    legacy_evaporation = vec(sum(legacy_soil.water.evaporation; dims = 1))
    aligned_evaporation = vec(sum(default_soil.water.evaporation; dims = 1))
    @test aligned_evaporation[1] ≈ 2.0f0 * legacy_evaporation[1] rtol = 2.0f-6
    @test aligned_evaporation[2] ≈ 2.0f0 * legacy_evaporation[2] rtol = 2.0f-6

    available_liquid_water = sum(
        base_soil.water.relative_content[1:2, 3] .*
        base_soil.water.holding_capacity_storage[1:2, 3] .-
        base_crop.fluxes.water.transpiration_layer[1:2, 3],
    )
    @test aligned_evaporation[3] ≈ available_liquid_water atol = 1.0f-7
    @test aligned_evaporation[3] <= available_liquid_water
    @test legacy_evaporation[3] > available_liquid_water

    ratio_without_ice = Agrocosm.compute_lpjml_managed_soil_evaporation_ratio(
        0.1f0, 1.0f0, available_liquid_water, available_liquid_water,
        200.0f0, 0.0f0,
    )
    ratio_with_ice = Agrocosm.compute_lpjml_managed_soil_evaporation_ratio(
        0.1f0, 1.0f0, available_liquid_water + 50.0f0,
        available_liquid_water, 200.0f0, 0.0f0,
    )
    @test ratio_with_ice > ratio_without_ice
    @test ratio_with_ice <= 1.0f0
end
