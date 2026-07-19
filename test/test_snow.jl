using Agrocosm
using Test

@testset "LPJ-compatible snow fluxes" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    weather = Agrocosm.init_weather(1, identity)

    # Cold precipitation accumulates as snow, followed by LPJmL's fixed 0.1 mm sublimation.
    weather.temp .= -5.0f0
    weather.prec .= 3.0f0
    snow!(soil, weather)

    @test weather.prec[1] == 0.0f0
    @test soil.snowpack[1] ≈ 2.9f0 atol = 1.0f-6
    @test soil.snow_sublimation[1] == 0.1f0
    @test soil.snowmelt[1] == 0.0f0
    @test soil.snow_runoff[1] == 0.0f0

    # Melt water is added to liquid precipitation before interception and infiltration.
    snow_before = soil.snowpack[1]
    precipitation = 2.0f0
    weather.temp .= 5.0f0
    weather.prec .= precipitation
    snow!(soil, weather)

    @test soil.snowmelt[1] > 0.0f0
    @test weather.prec[1] ≈ precipitation + soil.snowmelt[1] atol = 1.0f-6
    @test snow_before + precipitation ≈
          soil.snowpack[1] + weather.prec[1] + soil.snow_sublimation[1] atol = 1.0f-5

    # Snow above the configured capacity is exported as snow runoff.
    small_snowpack_params = LPJmLParams{Float32}(; maxsnowpack = 1.0f0)
    soil.snowpack .= 0.0f0
    weather.temp .= -5.0f0
    weather.prec .= 2.0f0
    snow!(soil, weather; lpjmlparams = small_snowpack_params)

    @test soil.snow_runoff[1] == 1.0f0
    @test soil.snowpack[1] ≈ 0.9f0 atol = 1.0f-6
    @test soil.snow_sublimation[1] == 0.1f0
    @test weather.prec[1] == 0.0f0
end
