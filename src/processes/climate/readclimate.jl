"""
readclimate!(climate, day, dailyWeather)

Read one day of climate forcing into runtime weather buffers.
"""
function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      day::Integer
)

    dailyWeather.temp = climate.temp[day, :]
    dailyWeather.prec = climate.prec[day, :]
    dailyWeather.swr = climate.sw[day, :]
    dailyWeather.lwr = climate.lw[day, :]
    if hasproperty(climate, :wind)
        dailyWeather.wind = climate.wind[day, :]
    else
        fill!(dailyWeather.wind, lpjmlparams.volatil_wind)
    end
    dailyWeather.annual_co2 = ppm2Pa(climate.co2[[div(day-1, 365) + 1]])

end

function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      CO2::AbstractArray{T},
                      day::Integer
) where {T <: AbstractFloat}

    dailyWeather.temp = climate.temp[day, :]
    dailyWeather.prec = climate.prec[day, :]
    dailyWeather.swr = climate.sw[day, :]
    dailyWeather.lwr = climate.lw[day, :]
    if hasproperty(climate, :wind)
        dailyWeather.wind = climate.wind[day, :]
    else
        fill!(dailyWeather.wind, lpjmlparams.volatil_wind)
    end
    dailyWeather.daily_co2 = ppm2Pa(CO2[day, :])

end
