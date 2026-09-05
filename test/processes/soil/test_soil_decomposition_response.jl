using Agrocosm
using Test

@testset "LPJmL soil decomposition response and rates" begin
    @test lpjmlparams.k_soil10.fast ≈ 0.04f0 / 365.0f0
    @test lpjmlparams.k_soil10.slow ≈ 0.001f0 / 365.0f0
    @test lpjmlparams.atmfrac == 0.5f0

    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.thermal.temperature .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.wilting_storage .= 0.0f0
    soil.water.free_water .= 0.0f0
    soil.water.relative_content[:, 1] .= 0.5f0
    soil.water.relative_content[:, 2] .= 0.25f0
    soil.water.available_ice_storage[:, 2] .= 50.0f0
    soil.surface_litter.temperature .= -20.0f0
    soil.surface_litter.water_capacity .= 1.0f0
    soil.surface_litter.water_storage .= 0.5f0

    soil_decomp_response!(state)

    # Cell 2 has half of its pore volume occupied by ice. Its 25 mm liquid
    # therefore has the same 50% liquid-pore saturation as cell 1's 50 mm.
    @test soil.decomposition.response[:, 1] ≈
          soil.decomposition.response[:, 2] rtol = 2.0f-6
    @test all(iszero, soil.decomposition.litter_response[1, :])
    @test soil.decomposition.litter_response[2, :] ≈
          soil.decomposition.response[1, :]
    @test soil.decomposition.litter_response[3, :] ≈
          soil.decomposition.response[1, :]

    crop = init_crop(1, identity)
    decay_soil = init_soil(1, soilparams.soildepth, identity)
    decay_state = test_model_state(crop, decay_soil)
    decay_soil.thermal.temperature .= 10.0f0
    decay_soil.water.saturation_storage .= 100.0f0
    decay_soil.water.holding_capacity_storage .= 100.0f0
    decay_soil.water.relative_content .= 0.5f0
    decay_soil.carbon.fast .= 100.0f0
    decay_soil.carbon.slow .= 100.0f0
    soil_carbon!(decay_state, decay_state)
    response = decay_soil.decomposition.response[1, 1]
    @test decay_soil.carbon.decomposed_fast[1, 1] ≈
        -100.0f0 * expm1(-0.04f0 / 365.0f0 * response) rtol = 2.0f-5
    @test decay_soil.carbon.decomposed_slow[1, 1] ≈
        -100.0f0 * expm1(-0.001f0 / 365.0f0 * response) rtol = 2.0f-4
end

@testset "Surface litter response is not capped with mineral soil" begin
    for T in (Float32, Float64)
        soil = init_soil(T, 3, T.(soilparams.soildepth), identity)
        state = test_model_state(soil)
        soil.thermal.temperature .= reshape(T[25, 10, -20], 1, 3)
        soil.surface_litter.temperature .= T[25, 10, 25]
        soil.surface_litter.water_capacity .= one(T)
        soil.surface_litter.water_storage .= T(0.6)
        soil.water.saturation_storage .= T(100)
        soil.water.holding_capacity_storage .= T(100)
        soil.water.wilting_storage .= zero(T)
        soil.water.free_water .= zero(T)
        soil.water.relative_content .= T(0.6)
        soil.carbon.litter .= T(100)
        soil.nitrogen.litter .= soil.carbon.litter ./ T(15)
        soil.carbon.litter_response .= T[0.97 / 365, 0.97 / 365, 0.3 / 365]

        soil_cn_decomposition!(state; lpjmlparams = LPJmLParams{T}())

        # Original LPJmL 6.1.9 NO_METHANE C oracle, also uncapped in 5.10.0.
        expected_response = T[2.9683781149847515, 0.9274415597600001, 0]
        @test soil.decomposition.litter_response[1, :] ≈ expected_response rtol = 2e-5
        @test soil.decomposition.response[1, :] ≈ min.(expected_response, one(T)) rtol = 2e-5
        @test soil.decomposition.litter_response[2, :] == soil.decomposition.response[1, :]
        @test soil.decomposition.litter_response[3, :] == soil.decomposition.response[1, :]
        @test soil.carbon.decomposed_litter[1, :] ≈
              T[0.7857533412465187, 0.24616728086105866, 0] rtol = 2e-5
        @test soil.nitrogen.decomposed_litter ≈ soil.carbon.decomposed_litter ./ T(15) rtol = 2e-5
    end
end
