using Agrocosm
using Test

@testset "LPJmL-style nitrification" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    Agrocosm.initialize_soil_class_parameters!(
        soil.properties, Int32[1, 9],
    )
    soil.properties.ph .= 6.5f0
    soil.nitrogen.ammonium[1, :] .= 0.8f0
    soil.water.relative_content[1, :] .= 0.6f0
    soil.water.holding_capacity_storage[1, :] .= 100.0f0
    soil.water.saturation_storage[1, :] .= 100.0f0
    soil.thermal.temperature[1, :] .= 20.0f0
    ammonium_before = copy(soil.nitrogen.ammonium[1, :])
    nitrate_before = copy(soil.nitrogen.nitrate[1, :])

    Agrocosm.launch_1D!(
        Agrocosm.nitrify_kernel!,
        soil.properties.ph,
        soil.properties.nitrification_a,
        soil.properties.nitrification_b,
        soil.properties.nitrification_c,
        soil.properties.nitrification_d,
        soil.nitrogen.ammonium,
        soil.nitrogen.nitrate,
        soil.water.relative_content,
        soil.water.holding_capacity_storage,
        soil.water.wilting_storage,
        soil.water.wilting_ice_fraction,
        soil.water.free_water,
        soil.water.saturation_storage,
        soil.thermal.temperature,
        soil.nitrogen.nitrification,
        soil.nitrogen.n2o_nitrification,
        (; lpjmlparams, soil_layers = 5),
    )

    gross = soil.nitrogen.nitrification[1, :]
    n2o = soil.nitrogen.n2o_nitrification[1, :]
    @test all(gross .> 0.0f0)
    @test all(gross .<= ammonium_before)
    @test gross[2] > gross[1]
    @test n2o ≈ lpjmlparams.k_2 * gross atol = 1.0f-7
    @test ammonium_before - soil.nitrogen.ammonium[1, :] ≈ gross atol = 1.0f-7
    @test soil.nitrogen.nitrate[1, :] - nitrate_before ≈ gross - n2o atol = 1.0f-7
end
