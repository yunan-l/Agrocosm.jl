using Agrocosm
using Test

@testset "LPJmL-style ammonia volatilization" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.properties.ph .= 7.0f0
    soil.nitrogen.ammonium[1, 1] = 1.0f0
    air_temperature = Float32[20.0]
    ammonium_before = soil.nitrogen.ammonium[1, 1]

    Agrocosm.launch_1D!(
        Agrocosm.volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        air_temperature,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )

    flux = soil.nitrogen.volatilization[1]
    @test 0.0f0 < flux <= ammonium_before
    @test ammonium_before - soil.nitrogen.ammonium[1, 1] ≈ flux atol = 1.0f-7
    @test all(soil.nitrogen.ammonium .>= 0.0f0)
end
