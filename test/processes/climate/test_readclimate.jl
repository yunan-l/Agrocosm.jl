using Agrocosm
using Test

@testset "Optional daily wind forcing" begin
    weather = init_weather(2, identity)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)

    climate_without_wind = (
        temp = Float32[10 11; 12 13],
        prec = zeros(Float32, 2, 2),
        sw = zeros(Float32, 2, 2),
        lw = zeros(Float32, 2, 2),
        co2 = Float32[400],
    )
    readclimate!(climate_without_wind, weather, 1)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)

    climate_with_wind = merge(
        climate_without_wind,
        (wind = Float32[2.0 3.0; 4.0 5.0],),
    )
    readclimate!(climate_with_wind, weather, 2)
    @test weather.wind == Float32[4, 5]

    # Moving back to a legacy archive must not retain the previous day's wind.
    readclimate!(climate_without_wind, weather, 2)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)

    loader_input = (
        temp_spinup = zeros(Float32, 2, 2),
        temp = zeros(Float32, 2, 2),
        prec = zeros(Float32, 2, 2),
        swdown = zeros(Float32, 2, 2),
        lwnet = zeros(Float32, 2, 2),
        co2 = Float32[400],
        temp_n = zeros(Float32, 2, 2),
        prec_n = zeros(Float32, 2, 2),
        sw_n = zeros(Float32, 2, 2),
        lw_n = zeros(Float32, 2, 2),
        windspeed = Float32[1 2; 3 4],
    )
    loaded = ClimateDataLoader(loader_input, [2], identity)
    @test hasproperty(loaded, :wind)
    @test loaded.wind == Float32[2; 4;;]
end
