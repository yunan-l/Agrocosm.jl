using Agrocosm
using Test

@testset "Daily freeze-thaw diagnostics" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    diagnostics = init_thermal_balance(1, 2, identity)
    soil.water.storage .= 50.0f0
    pedotransfer!(soil)
    soil_temperature!(soil, Float32[-20.0, 10.0], Float32[2.0, 2.0])
    Agrocosm.record_thermal_balance!(diagnostics, 1, soil)

    @test diagnostics.total_ice_storage[1, 1] > 0.0f0
    @test diagnostics.total_ice_storage[1, 2] == 0.0f0
    @test diagnostics.maximum_frozen_fraction[1, 1] > 0.0f0
    @test diagnostics.wilting_ice_storage[1, 1] > 0.0f0
    @test diagnostics.available_ice_storage[1, 1] > 0.0f0
    @test diagnostics.free_ice_storage[1, 1] >= 0.0f0
    @test diagnostics.ice_pool_residual[1, 1] <= 2.0f-5
    @test all(isfinite, diagnostics.energy_residual)
    @test diagnostics.minimum_temperature[1, 1] <=
          diagnostics.maximum_temperature[1, 1]
end
