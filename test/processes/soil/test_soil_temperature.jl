using Agrocosm
using Test

@testset "Soil-temperature baseline" begin
    # A perfect linear 31-day history must be recovered without an off-by-one.
    x = Float32.(1:31)
    climate_history = reshape(2.0f0 .+ 3.0f0 .* x, 31, 1)
    intercept, slope = Agrocosm.linreg(climate_history)
    @test intercept[1] ≈ 2.0f0 atol = 1.0f-4
    @test slope[1] ≈ 3.0f0 atol = 1.0f-4

    soil = init_soil(1, soilparams.soildepth, identity)
    climbuf = init_climbuf(1, identity)
    climbuf.temp[:, 1] .= Float32.(1:31)
    climbuf.atemp_mean .= 16.0f0
    soil.water.relative_content .= 0.5f0
    soil.thermal.diffusivity_0 .= 0.5f0
    soil.thermal.diffusivity_15 .= 0.6f0

    soiltemp_lag!(soil, climbuf)

    @test all(isfinite, soil.thermal.temperature)
    @test all(soil.thermal.temperature[:, 1] .== soil.thermal.temperature[1, 1])

    # The current dry-soil fallback assigns buffered air temperature to all layers.
    soil.water.relative_content[1, 1] = 0.0f0
    soiltemp_lag!(soil, climbuf)
    @test all(soil.thermal.temperature[:, 1] .== climbuf.temp[30, 1])
end

@testset "LPJmL soil decomposition temperature bounds" begin
    temperatures = Float32[-30.0, -15.01, -15.0, 10.0, 40.0, 55.0]
    response = Agrocosm.temp_response(temperatures)

    @test response[1] == 0.0f0
    @test response[2] == 0.0f0
    @test response[3] > 0.0f0
    @test response[4] ≈ 1.0f0 atol = 1.0f-6
    @test response[5] == response[6]
    @test all(isfinite, response)
end
