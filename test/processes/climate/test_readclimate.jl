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
        windspeed = Float32[1 2; 3 4],
    )
    loaded = ClimateDataLoader(loader_input, [2], identity)
    @test hasproperty(loaded, :wind)
    @test loaded.wind == Float32[2; 4;;]

    loaded64 = ClimateDataLoader(loader_input, [2], identity; T = Float64)
    @test eltype(loaded64.temp) == Float64
    @test eltype(loaded64.co2) == Float64
    @test eltype(loaded64.wind) == Float64
end


@testset "Climate read kernel matches vector reference" begin
    cells = 5
    days = 3
    climate = (
        temp = reshape(Float32.(1:(days * cells)), days, cells),
        prec = reshape(Float32.(11:(10 + days * cells)), days, cells),
        sw = reshape(Float32.(101:(100 + days * cells)), days, cells),
        lw = reshape(Float32.(-20:(-21 + days * cells)), days, cells),
        wind = reshape(Float32.(range(1, 6; length = days * cells)), days, cells),
        co2 = Float32[400, 405],
    )
    reference = init_weather(cells, identity)
    kernel = init_weather(cells, identity)
    reference_co2 = Agrocosm.readclimate_reference!(climate, reference, 2)
    kernel_co2 = readclimate!(climate, kernel, 2)
    for field in (:temp, :prec, :swr, :lwr, :wind, :annual_co2)
        @test getproperty(kernel, field) ≈ getproperty(reference, field)
    end
    @test reference_co2 === reference.annual_co2
    @test kernel_co2 === kernel.annual_co2

    daily_co2 = reshape(Float32.(range(390, 430; length = days * cells)), days, cells)
    daily_climate = merge(climate, (co2 = daily_co2,))
    reference_co2 = Agrocosm.readclimate_reference!(daily_climate, reference, 3)
    kernel_co2 = readclimate!(daily_climate, kernel, 3)
    for field in (:temp, :prec, :swr, :lwr, :wind, :daily_co2)
        @test getproperty(kernel, field) ≈ getproperty(reference, field)
    end
    @test reference_co2 === reference.daily_co2
    @test kernel_co2 === kernel.daily_co2
end
