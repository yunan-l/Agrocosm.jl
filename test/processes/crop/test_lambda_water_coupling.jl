using Agrocosm
using Test

@testset "LPJ CO2 units and minimum conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= Float32[100, 150, 250, 400, 400]
    pedotransfer!(state)

    crop.state.phenology.is_growing .= 1
    crop.state.carbon.root .= 1000.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.auxiliary.canopy.canopy_wet .= 0.0f0
    pet.eeq .= 1.0f0
    pet.daylength .= 12.0f0

    adtmm = Float32[5.0]
    co2_pa = Float32[40.0]
    expected_gp = 1.6f0 * adtmm[1] /
                  (co2_pa[1] * 1.0f-5 * (1.0f0 - 0.8f0) * 12.0f0 * 3600.0f0) +
                  cft1.gmin * crop.auxiliary.canopy.fpar[1]

    transpiration!(adtmm, cft1, state, pet, state, co2_pa)

    @test crop.auxiliary.canopy.canopy_conductance[1] ≈ expected_gp rtol = 1.0f-6
    @test crop.auxiliary.canopy.canopy_conductance[1] > cft1.gmin * crop.auxiliary.canopy.fpar[1]
end

@testset "Pre-phenology raw conductance drives the aligned water pass" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    crop.state.phenology.is_growing .= 1
    crop.state.carbon.root .= 1000.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    crop.auxiliary.canopy.fpar .= 0.2f0
    crop.fluxes.carbon.water_limited_assimilation .= 1.0f0
    pet.daylength .= 12.0f0
    co2 = Float32[40.0]

    Agrocosm.prepare_prephenology_canopy_conductance!(
        cft1, state, pet.daylength, co2,
    )
    raw_conductance = Agrocosm.compute_canopy_conductance(
        1.0f0, 40.0f0, 12.0f0, 0.2f0, cft1.gmin, lpjmlparams.LAMBDA_OPT,
    )
    @test crop.auxiliary.canopy.canopy_conductance[1] == raw_conductance

    # Mimic the post-phenology APAR/photosynthesis refresh. The aligned water
    # pass must retain raw `gp`; no current-fPAR/phenology multiplier is applied.
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.fluxes.carbon.water_limited_assimilation .= 10.0f0
    crop.auxiliary.canopy.canopy_wet .= 0.0f0
    soil.water.relative_content .= 0.0f0
    soil.water.relative_content[1] = 0.1f0
    soil.water.holding_capacity_storage .= 100.0f0
    pet.eeq .= 10.0f0

    recomputed_conductance = Agrocosm.compute_canopy_conductance(
        10.0f0, 40.0f0, 12.0f0, 0.8f0, cft1.gmin, lpjmlparams.LAMBDA_OPT,
    )
    @test recomputed_conductance != raw_conductance

    transpiration!(
        crop.fluxes.carbon.water_limited_assimilation,
        cft1, state, pet, state, co2;
        use_precomputed_conductance = true,
    )

    expected_demand = Agrocosm.compute_transpiration_demand(
        0.0f0, 10.0f0, lpjmlparams.ALPHAM, lpjmlparams.GM, raw_conductance,
    )
    expected_wscal = min(
        1.0f0,
        cft1.emax * 0.1f0 / expected_demand,
    )
    @test crop.state.water.demand_sum[1] ≈ expected_demand rtol = 1.0f-6
    @test crop.state.water.sufficiency[1] ≈ expected_wscal rtol = 1.0f-6
end

@testset "Inactive crop has zero conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    crop.state.phenology.is_growing .= 0
    crop.auxiliary.canopy.fpar .= 0.8f0
    pet.daylength .= 12.0f0
    pet.eeq .= 1.0f0

    transpiration!(Float32[5.0], cft1, state, pet, state, Float32[40.0])

    @test crop.auxiliary.canopy.canopy_conductance[1] == 0.0f0
end

@testset "Nitrogen limitation triggers only the LPJ water correction" begin
    @test Agrocosm.nitrogen_water_recoupling_required(
        2.0f0, 1.8f0, 3.0f0, 2.0f0,
    )
    @test !Agrocosm.nitrogen_water_recoupling_required(
        2.0f0, 1.995f0, 3.0f0, 2.0f0,
    )
    @test !Agrocosm.nitrogen_water_recoupling_required(
        2.0f0, 1.8f0, 2.05f0, 2.0f0,
    )
end

@testset "Nitrogen recoupling stores potential conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)
    crop.state.phenology.is_growing .= 1
    crop.auxiliary.photosynthesis.lambda .= 0.7f0
    crop.auxiliary.photosynthesis.temperature_stress .= 1.0f0
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.auxiliary.canopy.canopy_conductance .= 0.1f0
    crop.fluxes.carbon.water_limited_assimilation .= 6.0f0

    Agrocosm.refresh_potential_canopy_conductance!(
        cft1, state, Float32[12.0], Float32[40.0],
    )

    expected = Agrocosm.compute_canopy_conductance(
        6.0f0, 40.0f0, 12.0f0, 0.8f0, cft1.gmin, 0.7f0,
    )
    @test crop.auxiliary.canopy.canopy_conductance[1] ≈ expected
end

@testset "Nitrogen-limited transpiration uses final assimilation" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    crop.state.phenology.is_growing .= 1
    crop.auxiliary.photosynthesis.lambda .= 0.5f0
    crop.auxiliary.photosynthesis.temperature_stress .= 1.0f0
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.auxiliary.canopy.canopy_wet .= 0.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    crop.fluxes.carbon.water_limited_assimilation .= 5.0f0
    soil.water.relative_content .= 0.0f0
    soil.water.relative_content[1] = 0.5f0
    soil.water.holding_capacity_storage .= 100.0f0
    pet.daylength .= 12.0f0
    pet.eeq .= 2.0f0
    co2 = Float32[40.0]

    Agrocosm.finalize_nitrogen_limited_transpiration!(
        cft1, state, pet, state, co2,
    )

    expected_conductance = Agrocosm.compute_canopy_conductance(
        5.0f0, 40.0f0, 12.0f0, 0.8f0, cft1.gmin, 0.5f0,
    )
    expected_demand = Agrocosm.compute_transpiration_demand(
        0.0f0, 2.0f0, lpjmlparams.ALPHAM, lpjmlparams.GM,
        expected_conductance,
    )
    @test crop.auxiliary.canopy.canopy_conductance[1] ≈ expected_conductance
    @test crop.fluxes.water.transpiration_layer[1, 1] ≈ expected_demand * 0.8f0
    @test all(crop.fluxes.water.transpiration_layer[2:end, 1] .== 0.0f0)
end

@testset "Nitrogen-water recoupling resolves lambda within the old bracket" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    crop.state.phenology.is_growing .= 1
    crop.state.carbon.root .= 100.0f0
    crop.auxiliary.photosynthesis.lambda .= 0.7f0
    crop.auxiliary.photosynthesis.vcmax .= 1.0f0
    crop.auxiliary.photosynthesis.temperature_stress .= 1.0f0
    crop.auxiliary.canopy.canopy_conductance .= 10.0f0
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.auxiliary.canopy.apar .= 1.0f6
    crop.auxiliary.canopy.canopy_wet .= 0.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    crop.fluxes.carbon.water_limited_assimilation .= 6.0f0
    soil.water.relative_content .= 0.0f0
    soil.water.relative_content[1] = 0.8f0
    pet.daylength .= 12.0f0
    pet.eeq .= 20.0f0

    Agrocosm.recouple_nitrogen_water!(
        Val(:C3), cft1, state, pet, state, Float32[20.0], Float32[40.0],
    )

    @test 0.02f0 <= crop.auxiliary.photosynthesis.lambda[1] <= 0.7f0
    @test crop.auxiliary.photosynthesis.lambda[1] != 0.7f0

    crop.auxiliary.photosynthesis.lambda .= 0.7f0
    crop.auxiliary.canopy.canopy_conductance .= 10.0f0
    Agrocosm.recouple_nitrogen_water!(
        Val(:C4), cft3, state, pet, state, Float32[25.0], Float32[40.0],
    )
    @test 0.02f0 <= crop.auxiliary.photosynthesis.lambda[1] <= 0.7f0
end
