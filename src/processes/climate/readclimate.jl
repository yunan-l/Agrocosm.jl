"""
readclimate!(climate, day, dailyWeather)

Read one day of climate forcing into runtime weather buffers.
"""
function readclimate_reference!(climate::NamedTuple,
                                dailyWeather::DailyWeather,
                                day::Integer
)

    @views dailyWeather.temp .= climate.temp[day, :]
    @views dailyWeather.prec .= climate.prec[day, :]
    @views dailyWeather.swr .= climate.sw[day, :]
    @views dailyWeather.lwr .= climate.lw[day, :]
    if hasproperty(climate, :wind)
        @views dailyWeather.wind .= climate.wind[day, :]
    else
        fill!(dailyWeather.wind, lpjmlparams.volatil_wind)
    end
    co2_year = div(day - 1, 365) + 1
    @views dailyWeather.annual_co2 .= climate.co2[co2_year:co2_year] .* 0.1f0

end

function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      day::Integer)
    has_wind = hasproperty(climate, :wind)
    wind = has_wind ? climate.wind : climate.temp
    launch_1D!(
        read_annual_climate_kernel!,
        dailyWeather.temp,
        dailyWeather.prec,
        dailyWeather.swr,
        dailyWeather.lwr,
        dailyWeather.wind,
        dailyWeather.annual_co2,
        climate.temp,
        climate.prec,
        climate.sw,
        climate.lw,
        wind,
        climate.co2,
        day,
        has_wind,
        eltype(dailyWeather.temp)(lpjmlparams.volatil_wind),
    )
    return nothing
end

function readclimate_reference!(climate::NamedTuple,
                                dailyWeather::DailyWeather,
                                CO2::AbstractArray{T},
                                day::Integer
) where {T <: AbstractFloat}

    @views dailyWeather.temp .= climate.temp[day, :]
    @views dailyWeather.prec .= climate.prec[day, :]
    @views dailyWeather.swr .= climate.sw[day, :]
    @views dailyWeather.lwr .= climate.lw[day, :]
    if hasproperty(climate, :wind)
        @views dailyWeather.wind .= climate.wind[day, :]
    else
        fill!(dailyWeather.wind, lpjmlparams.volatil_wind)
    end
    @views dailyWeather.daily_co2 .= CO2[day, :] .* T(0.1)

end


function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      CO2::AbstractArray{T},
                      day::Integer) where {T <: AbstractFloat}
    has_wind = hasproperty(climate, :wind)
    wind = has_wind ? climate.wind : climate.temp
    launch_1D!(
        read_daily_climate_kernel!,
        dailyWeather.temp,
        dailyWeather.prec,
        dailyWeather.swr,
        dailyWeather.lwr,
        dailyWeather.wind,
        dailyWeather.daily_co2,
        climate.temp,
        climate.prec,
        climate.sw,
        climate.lw,
        wind,
        CO2,
        day,
        has_wind,
        T(lpjmlparams.volatil_wind),
    )
    return nothing
end

@kernel inbounds = true function read_annual_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    annual_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    co2_forcing::AbstractVector{T},
    day::Integer,
    has_wind::Bool,
    default_wind::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    if cell == 1
        co2_year = div(day - 1, 365) + 1
        annual_co2[1] = co2_forcing[co2_year] * T(0.1)
    end
end

@kernel inbounds = true function read_daily_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    daily_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    co2_forcing::AbstractMatrix{T},
    day::Integer,
    has_wind::Bool,
    default_wind::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    daily_co2[cell] = co2_forcing[day, cell] * T(0.1)
end
