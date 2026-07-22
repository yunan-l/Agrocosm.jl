"""
snow!(soil, dailyWeather)

Update snowpack, snow height, and snow cover fraction from daily temperature
and precipitation forcing.
"""
function snow!(soil,
               dailyWeather::DailyWeather;
               snowparams::SnowParams = snowparams,
               lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (; snowparams, lpjmlparams)

    launch_1D!(
        snow_kernel!,
        dailyWeather.temp,
        dailyWeather.prec,
        soil_snow_prognostic(soil).pack,
        soil_snow_fluxes(soil).melt,
        soil_snow_fluxes(soil).sublimation,
        soil_snow_fluxes(soil).runoff,
        soil_snow_prognostic(soil).height,
        soil_snow_prognostic(soil).fraction,
        kernel_params
    )

end

function snow!(state::ModelState;
               snowparams::SnowParams = snowparams,
               lpjmlparams::LPJmLParams = lpjmlparams)
    weather = state.inputs.weather
    snow_state = state.prognostic.soil.snow
    snow_fluxes = state.fluxes.soil.snow
    launch_1D!(
        snow_kernel!, weather.temp, weather.prec,
        snow_state.pack, snow_fluxes.melt, snow_fluxes.sublimation,
        snow_fluxes.runoff, snow_state.height, snow_state.fraction,
        (; snowparams, lpjmlparams),
    )
    return nothing
end



@kernel inbounds = true function snow_kernel!(
                              temp::AbstractArray{T},
                              prec::AbstractArray{T},
                              soil_snowpack::AbstractArray{T},
                              soil_snowmelt::AbstractArray{T},
                              soil_snow_sublimation::AbstractArray{T},
                              soil_snow_runoff::AbstractArray{T},
                              soil_snowheight::AbstractArray{T},
                              soil_snowfraction::AbstractArray{T},
                              kernel_params
) where {T <: AbstractFloat}

    cell = @index(Global)

    @unpack tsnow, snow_skin_depth, th_diff_snow, lambda_snow, c_water2ice, c_watertosnow, c_roughness= kernel_params.snowparams
    @unpack maxsnowpack = kernel_params.lpjmlparams

    soil_snowmelt[cell] = zero(T)
    soil_snow_sublimation[cell] = zero(T)
    soil_snow_runoff[cell] = zero(T)

    # precipitation falls as snow
    if temp[cell] < tsnow
        soil_snowpack[cell] += prec[cell]
        if soil_snowpack[cell] > maxsnowpack
            soil_snow_runoff[cell] = soil_snowpack[cell] - maxsnowpack
            soil_snowpack[cell] = maxsnowpack
        end
        prec[cell] = zero(T)
    end

    # sublimation of snow
    if soil_snowpack[cell] > T(0.1)
        soil_snowpack[cell] -= T(0.1)
        soil_snow_sublimation[cell] = T(0.1)
    end

    # snow layer is insulating
    timestep2sec = T(24.0 * 3600.0)
    if soil_snowpack[cell] > T(1e-7)
        if temp[cell] > zero(T)
            depth = min(soil_snowpack[cell], snow_skin_depth)
            dT = th_diff_snow * timestep2sec / (depth * depth) * T(1000000.0) * (temp[cell] - tsnow)
            heatflux = lambda_snow * (tsnow - zero(T) + dT) / depth * T(1000)
            melt_heat = min(heatflux * timestep2sec, depth * T(1e-3) * c_water2ice) #[J/m2]
            melt = melt_heat / c_water2ice * T(1000)
            soil_snowmelt[cell] += melt
            soil_snowpack[cell] -= melt
            if soil_snowpack[cell] < T(1e-7)
                soil_snowpack[cell] = zero(T)
                soil_snowheight[cell] = zero(T)
                soil_snowfraction[cell] = zero(T)
            end
        end
    end

    # Add melt water to rainfall before interception and infiltration.
    prec[cell] += soil_snowmelt[cell]

    # calculate snow height and fraction of snow coverage
    if soil_snowpack[cell] > T(1e-7)
        HS = c_watertosnow * (soil_snowpack[cell] / T(1000.0)) # mm -> m */
        frsg = HS / (HS+ T(0.5) * c_roughness)
        soil_snowheight[cell] = HS
        soil_snowfraction[cell] = frsg
    else
        soil_snowpack[cell] = zero(T)
        soil_snowheight[cell] = zero(T)
        soil_snowfraction[cell] = zero(T)
    end

end
