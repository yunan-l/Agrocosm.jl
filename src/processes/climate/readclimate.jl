"""
    readclimate!(climate, dailyWeather, day)

Read one day of climate forcing and return the active CO₂ buffer. A vector
`climate.co2` is interpreted as a global annual series shared by all cells; a
matrix is interpreted as daily forcing with shape `(day, cell)`.
"""

function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      day::Integer)
    has_wind = hasproperty(climate, :wind)
    has_no3_deposition = hasproperty(climate, :no3_deposition)
    has_nh4_deposition = hasproperty(climate, :nh4_deposition)
    co2_daily = hasproperty(climate, :co2_daily) && climate.co2_daily
    wind = has_wind ? climate.wind : climate.temp
    no3_deposition = has_no3_deposition ? climate.no3_deposition : climate.temp
    nh4_deposition = has_nh4_deposition ? climate.nh4_deposition : climate.temp
    default_wind = eltype(dailyWeather.temp)(lpjmlparams.volatil_wind)
    read_vapour_pressure!(climate, dailyWeather, day)
    if ndims(climate.co2) == 1
        launch_1D!(
            read_annual_climate_kernel!,
            dailyWeather.temp,
            dailyWeather.prec,
            dailyWeather.swr,
            dailyWeather.lwr,
            dailyWeather.wind,
            dailyWeather.no3_deposition,
            dailyWeather.nh4_deposition,
            dailyWeather.annual_co2,
            climate.temp,
            climate.prec,
            climate.sw,
            climate.lw,
            wind,
            no3_deposition,
            nh4_deposition,
            climate.co2,
            day,
            has_wind,
            has_no3_deposition,
            has_nh4_deposition,
            co2_daily,
            default_wind,
        )
        # AFTER the precipitation is filled, not before: called earlier this
        # carried yesterday's rain into today's interception, which showed up as
        # a one-bit difference in the checkpoint test and nowhere else.
        read_canopy_rain!(climate, dailyWeather, day)
        return dailyWeather.annual_co2
    elseif ndims(climate.co2) == 2
        launch_1D!(
            read_daily_climate_kernel!,
            dailyWeather.temp,
            dailyWeather.prec,
            dailyWeather.swr,
            dailyWeather.lwr,
            dailyWeather.wind,
            dailyWeather.no3_deposition,
            dailyWeather.nh4_deposition,
            dailyWeather.daily_co2,
            climate.temp,
            climate.prec,
            climate.sw,
            climate.lw,
            wind,
            no3_deposition,
            nh4_deposition,
            climate.co2,
            day,
            has_wind,
            has_no3_deposition,
            has_nh4_deposition,
            default_wind,
        )
        read_canopy_rain!(climate, dailyWeather, day)
        return dailyWeather.daily_co2
    else
        throw(ArgumentError("climate.co2 must be a vector or a (day, cell) matrix"))
    end
end

"""
    read_canopy_rain!(climate, dailyWeather, day)

Set the water that reaches the CANOPY: the day's precipitation, less the
irrigation when the run says its irrigation goes below the canopy.

Basin and furrow irrigation do not wet a canopy. Entering them as rainfall does,
and `canopy_wet` is proportional to the day's rain while transpiration demand
carries a factor `1 - canopy_wet`, so a 317 mm basin irrigation takes that day's
transpiration to nearly zero. Measured at Maricopa: the WET arm's mean demand
comes out BELOW the dry arm's, 2.739 against 2.867 mm/day, because it is
irrigated more often.
"""
function read_canopy_rain!(climate::NamedTuple,
                           dailyWeather::DailyWeather,
                           day::Integer)
    subcanopy = hasproperty(climate, :subcanopy_irrigation) &&
                climate.subcanopy_irrigation && hasproperty(climate, :irrigation)
    if !subcanopy
        copyto!(dailyWeather.canopy_rain, dailyWeather.prec)
        return nothing
    end
    launch_1D!(
        read_canopy_rain_kernel!,
        dailyWeather.canopy_rain, dailyWeather.prec, climate.irrigation, day,
    )
    return nothing
end

@kernel inbounds = true function read_canopy_rain_kernel!(
    canopy_rain::AbstractVector{T},
    precipitation::AbstractVector{T},
    irrigation::AbstractMatrix{T},
    day::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    canopy_rain[cell] = max(precipitation[cell] - irrigation[day, cell], zero(T))
end

"""
    read_vapour_pressure!(climate, dailyWeather, day)

Fill today's actual vapour pressure (kPa) from specific humidity and surface
pressure, or leave it at zero when the forcing carries neither.

Zero is the signal that an aerodynamically coupled transpiration demand cannot
be formed, and the process falls back to LPJmL's uncoupled one rather than
inventing a humidity. GSWP3-W5E5 ships `huss`, so this is a wiring question and
not a data one.
"""
function read_vapour_pressure!(climate::NamedTuple,
                               dailyWeather::DailyWeather,
                               day::Integer)
    if !(hasproperty(climate, :specific_humidity) &&
         hasproperty(climate, :surface_pressure) &&
         hasproperty(climate, :diurnal_range))
        fill!(dailyWeather.vapour_deficit, zero(eltype(dailyWeather.vapour_deficit)))
        return nothing
    end
    launch_1D!(
        read_vapour_pressure_kernel!,
        dailyWeather.vapour_deficit,
        climate.specific_humidity,
        climate.surface_pressure,
        climate.temp,
        climate.diurnal_range,
        day,
    )
    return nothing
end

@kernel inbounds = true function read_vapour_pressure_kernel!(
    vapour_deficit::AbstractVector{T},
    humidity_forcing::AbstractMatrix{T},
    pressure_forcing::AbstractMatrix{T},
    temperature_forcing::AbstractMatrix{T},
    diurnal_range_forcing::AbstractMatrix{T},
    day::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    humidity = humidity_forcing[day, cell]
    pressure = pressure_forcing[day, cell]
    actual = vapour_pressure_from_specific_humidity(humidity, pressure)
    mean_temperature = temperature_forcing[day, cell]
    half_range = max(diurnal_range_forcing[day, cell], zero(T)) / T(2)
    saturated = (saturation_vapour_pressure(mean_temperature + half_range) +
                 saturation_vapour_pressure(mean_temperature - half_range)) / T(2)
    vapour_deficit[cell] = max(saturated - actual, zero(T))
end

@kernel inbounds = true function read_annual_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    no3_deposition::AbstractVector{T},
    nh4_deposition::AbstractVector{T},
    annual_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    no3_deposition_forcing::AbstractMatrix{T},
    nh4_deposition_forcing::AbstractMatrix{T},
    co2_forcing::AbstractVector{T},
    day::Integer,
    has_wind::Bool,
    has_no3_deposition::Bool,
    has_nh4_deposition::Bool,
    co2_daily::Bool,
    default_wind::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    no3_deposition[cell] = has_no3_deposition ? no3_deposition_forcing[day, cell] : zero(T)
    nh4_deposition[cell] = has_nh4_deposition ? nh4_deposition_forcing[day, cell] : zero(T)
    if cell == 1
        co2_index = co2_daily ? day : div(day - 1, 365) + 1
        annual_co2[1] = co2_forcing[co2_index] * T(0.1)
    end
end

@kernel inbounds = true function read_daily_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    no3_deposition::AbstractVector{T},
    nh4_deposition::AbstractVector{T},
    daily_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    no3_deposition_forcing::AbstractMatrix{T},
    nh4_deposition_forcing::AbstractMatrix{T},
    co2_forcing::AbstractMatrix{T},
    day::Integer,
    has_wind::Bool,
    has_no3_deposition::Bool,
    has_nh4_deposition::Bool,
    default_wind::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    no3_deposition[cell] = has_no3_deposition ? no3_deposition_forcing[day, cell] : zero(T)
    nh4_deposition[cell] = has_nh4_deposition ? nh4_deposition_forcing[day, cell] : zero(T)
    daily_co2[cell] = co2_forcing[day, cell] * T(0.1)
end
