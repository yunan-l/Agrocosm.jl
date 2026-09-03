using Agrocosm
using Test

@testset "LPJmL-style ammonia volatilization" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.properties.ph .= 7.0f0
    soil.nitrogen.ammonium[1, 1] = 1.0f0
    air_temperature = Float32[20.0]
    # Use low, non-saturating wind speeds so the monotonic wind response is
    # observable instead of both cases being capped by the available NH4.
    wind_speed = Float32[0.01]
    ammonium_before = soil.nitrogen.ammonium[1, 1]

    Agrocosm.launch_1D!(
        Agrocosm.volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        air_temperature,
        wind_speed,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )

    flux = soil.nitrogen.volatilization[1]
    @test 0.0f0 < flux <= ammonium_before
    @test ammonium_before - soil.nitrogen.ammonium[1, 1] ≈ flux atol = 1.0f-7
    @test all(soil.nitrogen.ammonium .>= 0.0f0)

    soil.nitrogen.ammonium[1, 1] = ammonium_before
    high_wind = Float32[0.02]
    Agrocosm.launch_1D!(
        Agrocosm.volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        air_temperature,
        high_wind,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )
    @test soil.nitrogen.volatilization[1] > flux
end

@testset "LPJmL Float32 ammonia dissociation" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.properties.ph .= 6.6f0
    soil.nitrogen.ammonium[1, 1] = 9.8688955f0
    air_temperature = Float32[20.9]
    wind_speed = Float32[1.5]

    post_crop_nitrogen_losses!(
        test_model_state(soil);
        air_temperature,
        wind_speed,
        exact_lpjml_volatilization = true,
    )

    ammonium = 9.8688955
    temperature = 20.9
    ph = 6.6
    kelvin = temperature + 273.15
    dissociation = 10.0^(0.05 - 2788.0 / kelvin)
    aqueous_fraction = 1.0 / (1.0 + 10.0^(-ph) / dissociation)
    aqueous_nh3 = aqueous_fraction * ammonium / 200.0 * 1000.0
    henry = 0.2138 / kelvin * 10.0^(6.123 - 1825.0 / kelvin)
    transfer = 0.000612 * 1.5^0.8 * kelvin^0.382
    expected = Float32(86400.0 * transfer * henry * aqueous_nh3)

    @test soil.nitrogen.volatilization[1] ≈ expected rtol = 5.0f-6
    @test soil.nitrogen.ammonium[1, 1] ≈ 9.8688955f0 - expected rtol = 5.0f-6
end
