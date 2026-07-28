_pathway_albedo!(::Val{:C3}, pft, crop, soil, pet, maize) =
    albedo!(pft, crop, soil, pet)
_pathway_albedo!(::Val{:C4}, pft, crop, soil, pet, maize) =
    albedo!(pft, crop, soil, pet; maize)

_pathway_apar!(::Val{:C3}, pft, crop, pet, snow_height, maize) =
    apar_crop!(pft, crop, pet, snow_height)
_pathway_apar!(::Val{:C4}, pft, crop, pet, snow_height, maize) =
    apar_crop!(pft, crop, pet, snow_height; maize)

"""Execute one daily crop pathway over the canonical lifecycle state."""
function _daily_crop!(
    pathway::Union{Val{:C3}, Val{:C4}}, start_day, end_day,
    processes::ProcessModules, climate, state::ModelState;
    maize = true,
    irrigation = false,
    manure = false,
    fertilizer = :auto,
    with_tillage = true,
    nitrogen_limit_vcmax = false,
    water_balance = nothing,
    nitrogen_balance = nothing,
    carbon_balance = nothing,
    thermal_balance = nothing,
    simulation_day_offset::Integer = 0,
    diagnostic_offset::Integer = 0,
    reuse_output::Bool = false,
    selected_output::Union{Nothing, Set{Tuple{Symbol, Symbol}}} = nothing,
    prescribed_phu = nothing,
    prescribed_winter_type = nothing,
)
    pftparameters = processes.crop
    model_parameters = processes.global_parameters
    climbuf = state.prognostic.climate
    crop = state
    pet = state.auxiliary.pet
    soil = state
    managed_land = state.inputs.management
    dailyWeather = state.inputs.weather
    output = state.output

    T = eltype(crop_prognostic(crop).canopy.lai)
    pftparameters = convert_precision(T, pftparameters)
    model_parameters = convert_precision(T, model_parameters)
    global_params = model_parameters.lpjml
    photo_params = model_parameters.photosynthesis
    snow_params = model_parameters.snow
    thermal_params = model_parameters.soil_thermal
    decomp_params = model_parameters.soil_decomposition
    fertilizer = fertilizer_mode(fertilizer)
    automatic_fertilizer = fertilizer === :auto

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    annual_rows = count(
        climate_day -> (climate_day + simulation_day_offset) % 365 == 0,
        start_day:end_day,
    )
    output_rows = prepare_output_block!(
        output, end_day - start_day + 1, annual_rows;
        reuse = reuse_output, selected = selected_output,
    )
    annual_output_offset = 0

    for climate_day in start_day:end_day
        day = climate_day + simulation_day_offset
        block_day = climate_day - start_day + 1
        diagnostic_day = diagnostic_offset + block_day
        output_row = output_rows.first_daily_row + block_day - 1
        day_of_year = day % 365 != 0 ? day % 365 : 365
        current_co2 = readclimate!(climate, dailyWeather, climate_day)

        if carbon_balance !== nothing
            record_carbon_balance_start!(carbon_balance, diagnostic_day, crop, soil)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_start!(nitrogen_balance, diagnostic_day, crop, soil)
        end
        if water_balance !== nothing
            record_water_balance_start!(
                water_balance, diagnostic_day, soil, dailyWeather.prec,
            )
        end

        # Today's climate must enter history before sowing decisions.
        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day)
        cultivate!(
            crop, managed_land, soil, day_of_year;
            manure,
            apply_prescribed_fertilizer = fertilizer === :yes,
            prescribed_phu,
            prescribed_winter_type,
            lpjmlparams = global_params,
            laimax = pftparameters.laimax,
        )
        if carbon_balance !== nothing
            record_carbon_balance_after_cultivate!(carbon_balance, diagnostic_day, crop)
        end

        if with_tillage
            litter_tillage!(soil, crop)
            tillage_hydraulics!(soil, crop; lpjmlparams = global_params)
        end
        litter_bioturbation!(soil; lpjmlparams = global_params)

        # Radiation uses the snow state present at the start of the day.
        _pathway_albedo!(pathway, pftparameters, crop, soil, pet, maize)
        petpar!(
            pet, day_of_year, managed_land.latitude, dailyWeather.temp,
            dailyWeather.lwr, dailyWeather.swr,
        )
        snow!(soil, dailyWeather; snowparams = snow_params, lpjmlparams = global_params)
        if water_balance !== nothing
            record_water_balance_after_snow!(
                water_balance, diagnostic_day, dailyWeather.prec,
            )
        end

        pedotransfer!(soil; lpjmlparams = global_params)
        update_surface_litter_properties!(soil; thermalparams = thermal_params)
        soil_temperature!(
            soil, dailyWeather.temp, climbuf.atemp_mean;
            thermalparams = thermal_params, snowparams = snow_params,
        )

        # Existing litter/SOM decomposes before today's crop uptake.
        soil_cn_decomposition!(
            soil;
            lpjmlparams = global_params,
            soil_decomp_params = decomp_params,
        )

        phenology_crop!(
            crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength,
        )
        annual_output_row = day_of_year == 365 ?
            output_rows.first_annual_row + annual_output_offset : nothing
        harvest_crop!(
            crop, soil, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(soil, crop)
        annual_output_offset += day_of_year == 365
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, crop, soil,
                managed_land.residue_fraction,
            )
        end

        interception!(
            crop, pftparameters, pet.eeq, dailyWeather.prec;
            lpjmlparams = global_params,
        )
        pedotransfer!(soil; lpjmlparams = global_params)
        soil_infiltration!(
            soil, crop, dailyWeather.prec;
            irrigation,
            snowmelt = soil_snow_fluxes(soil).melt,
            air_temperature = dailyWeather.temp,
            lpjmlparams = global_params,
            thermalparams = thermal_params,
        )
        if thermal_balance !== nothing
            record_thermal_balance!(thermal_balance, diagnostic_day, soil)
        end

        _pathway_apar!(
            pathway, pftparameters, crop, pet,
            soil_snow_prognostic(soil).height, maize,
        )
        temp_stress(
            pftparameters, pet, crop, dailyWeather.temp;
            photoparams = photo_params,
        )
        photosynthesis!(
            pathway, pftparameters, crop, crop_canopy_auxiliary(crop).apar,
            pet.daylength, dailyWeather.temp, current_co2;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )

        transpiration!(
            crop_fluxes(crop).carbon.water_limited_assimilation,
            pftparameters, crop, pet, soil, current_co2;
            lpjmlparams = global_params,
        )
        solve_lambda!(
            pathway, pftparameters, crop, pet, dailyWeather.temp, current_co2;
            lpjmlparams = global_params,
            photoparams = photo_params,
        )

        if nitrogen_limit_vcmax
            crop_nitrogen!(
                crop, pftparameters, soil,
                crop_photosynthesis_auxiliary(crop).potential_vcmax,
                dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                lpjmlparams = global_params,
            )
            limit_vcmax_by_nitrogen!(
                crop, pftparameters, dailyWeather.temp;
                lpjmlparams = global_params,
            )
        end
        photosynthesis!(
            pathway, pftparameters, crop, crop_canopy_auxiliary(crop).apar,
            pet.daylength, dailyWeather.temp, current_co2;
            comp_vcmax = false,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )

        crop_carbon!(
            crop, output, pftparameters, dailyWeather.temp,
            soil_thermal_prognostic(soil).temperature;
            output_row, lpjmlparams = global_params,
        )
        terminate_failed_crop!(
            crop, soil, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(soil, crop)
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, crop, soil,
                managed_land.residue_fraction,
            )
        end

        if nitrogen_limit_vcmax
            allocate_crop_nitrogen!(crop, pftparameters)
        else
            crop_nitrogen!(
                crop, pftparameters, soil,
                crop_photosynthesis_auxiliary(crop).vcmax, dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                lpjmlparams = global_params,
            )
        end

        evaporation!(pet.eeq, crop, soil; lpjmlparams = global_params)
        soil_evapotranspiration!(soil, crop; irrigation)
        post_crop_nitrogen_losses!(
            soil;
            air_temperature = dailyWeather.temp,
            wind_speed = dailyWeather.wind,
            lpjmlparams = global_params,
        )

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, soil, crop)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_end!(nitrogen_balance, diagnostic_day, crop, soil)
        end
        if carbon_balance !== nothing
            record_carbon_balance_end!(carbon_balance, diagnostic_day, crop, soil)
        end
    end
    return nothing
end

"""Compatibility entry point for the C3-specialized daily crop transition."""
daily_crop_C3!(args...; kwargs...) = _daily_crop!(Val(:C3), args...; kwargs...)

"""Compatibility entry point for the C4-specialized daily crop transition."""
daily_crop_C4!(args...; kwargs...) = _daily_crop!(Val(:C4), args...; kwargs...)
