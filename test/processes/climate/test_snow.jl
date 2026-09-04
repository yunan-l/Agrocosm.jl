using Agrocosm
using Test

@testset "LPJ-compatible snow fluxes" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    weather = Agrocosm.init_weather(1, identity)
    state = test_model_state(soil; weather)

    # Cold precipitation accumulates as snow, followed by LPJmL's fixed 0.1 mm sublimation.
    weather.temp .= -5.0f0
    weather.prec .= 3.0f0
    snow!(state, weather)

    @test weather.prec[1] == 0.0f0
    @test soil.snow.pack[1] ≈ 2.9f0 atol = 1.0f-6
    @test soil.snow.sublimation[1] == 0.1f0
    @test soil.snow.melt[1] == 0.0f0
    @test soil.snow.runoff[1] == 0.0f0

    # LPJmL keeps melt separate while liquid rain passes through the canopy.
    snow_before = soil.snow.pack[1]
    precipitation = 2.0f0
    weather.temp .= 5.0f0
    weather.prec .= precipitation
    snow!(state, weather)

    melt = soil.snow.melt[1]
    @test melt > 0.0f0
    @test weather.prec[1] == precipitation
    @test snow_before + precipitation ≈
          soil.snow.pack[1] + weather.prec[1] + melt +
          soil.snow.sublimation[1] atol = 1.0f-5

    # Melt is restored only after canopy interception, before litter and soil.
    Agrocosm.add_snowmelt_to_precipitation!(weather.prec, soil.snow.melt)
    @test weather.prec[1] ≈ precipitation + melt atol = 1.0f-6

    # Snow above the configured capacity is exported as snow runoff.
    small_snowpack_params = LPJmLParams{Float32}(; maxsnowpack = 1.0f0)
    soil.snow.pack .= 0.0f0
    weather.temp .= -5.0f0
    weather.prec .= 2.0f0
    snow!(state, weather; lpjmlparams = small_snowpack_params)

    @test soil.snow.runoff[1] == 1.0f0
    @test soil.snow.pack[1] ≈ 0.9f0 atol = 1.0f-6
    @test soil.snow.sublimation[1] == 0.1f0
    @test weather.prec[1] == 0.0f0
end
