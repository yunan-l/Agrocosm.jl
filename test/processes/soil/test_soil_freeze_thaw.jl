using Agrocosm
using Test

@testset "LPJmL enthalpy-temperature phase relation" begin
    cf = 2.0f6
    cu = 3.0f6
    latent = 6.0f7
    @test Agrocosm.enthalpy_temperature(-4.0f6, cf, cu, latent) == -2.0f0
    @test Agrocosm.enthalpy_temperature(0.25f0 * latent, cf, cu, latent) == 0.0f0
    @test Agrocosm.enthalpy_temperature(latent + 6.0f6, cf, cu, latent) == 2.0f0
    @test Agrocosm.enthalpy_frozen_fraction(0.25f0 * latent, latent) == 0.75f0
end

@testset "Five-layer freeze-thaw conserves water and energy" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.water.storage .= 50.0f0
    pedotransfer!(soil)
    total_water_before = vec(sum(soil.water.storage + soil.water.ice_storage; dims = 1))

    soil_temperature!(soil, Float32[-20.0], Float32[2.0])
    total_water_after_freezing = vec(sum(
        soil.water.storage + soil.water.ice_storage; dims = 1,
    ))
    @test total_water_after_freezing ≈ total_water_before atol = 2.0f-5
    @test all((0.0f0 .<= soil.thermal.frozen_fraction) .&
              (soil.thermal.frozen_fraction .<= 1.0f0))
    @test maximum(soil.water.ice_storage) > 0.0f0
    phase_layers = (soil.thermal.frozen_fraction .> 0.0f0) .&
                   (soil.thermal.frozen_fraction .< 1.0f0)
    @test any(phase_layers)
    @test all(abs.(soil.thermal.temperature[phase_layers]) .<= 2.0f-5)
    relative_energy_residual = abs(soil.thermal.energy_residual[1]) /
        max(abs(soil.thermal.surface_energy_flux[1]), 1.0f0)
    @test relative_energy_residual < 2.0f-5

    ice_after_freezing = sum(soil.water.ice_storage)
    for _ in 1:5
        pedotransfer!(soil)
        soil_temperature!(soil, Float32[20.0])
    end
    @test sum(soil.water.ice_storage) < ice_after_freezing
    @test vec(sum(soil.water.storage + soil.water.ice_storage; dims = 1)) ≈
          total_water_before atol = 3.0f-5
    @test all(isfinite, soil.thermal.enthalpy)
end

@testset "Ice-liquid phase transfer closes water balance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    balance = init_water_balance(1, 1, identity)
    soil.water.storage .= 40.0f0
    pedotransfer!(soil)
    Agrocosm.record_water_balance_start!(balance, 1, soil, Float32[0.0])
    Agrocosm.record_water_balance_after_snow!(balance, 1, Float32[0.0])

    soil_temperature!(soil, Float32[-15.0], Float32[-2.0])
    Agrocosm.record_water_balance_end!(balance, 1, soil, crop)

    @test balance.soil_ice_storage_after[1, 1] > 0.0f0
    @test balance.residual[1, 1] ≈ 0.0f0 atol = 3.0f-5
end
