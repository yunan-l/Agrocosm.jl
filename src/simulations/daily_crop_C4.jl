"""Execute daily C4 processes over the canonical lifecycle state."""
function daily_crop_C4!(day_start, day_end,
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
        day_start:day_end,
    )
    output_rows = prepare_output_block!(output, day_end - day_start + 1, annual_rows)
    annual_output_offset = 0

    for climate_day = day_start : day_end

        day = climate_day + simulation_day_offset
        block_day = climate_day - day_start + 1
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
            record_water_balance_start!(water_balance, diagnostic_day, soil, dailyWeather.prec)
        end

        # Keep today's climate in history before sowing decisions, following
        # LPJmL's daily cell order and the dependency needed by dynamic sowing.
        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day)

        # initial crop variables in sowing day and fertilizer
        cultivate!(
            crop,
            managed_land,
            soil,
            day_of_year;
            manure = manure,
            apply_prescribed_fertilizer = fertilizer === :yes,
            lpjmlparams = global_params,
            laimax = pftparameters.laimax,
        )

        if carbon_balance !== nothing
            record_carbon_balance_after_cultivate!(carbon_balance, diagnostic_day, crop)
        end

        # LPJmL tills existing litter and loosens the topsoil at cultivation,
        # then applies daily agtop -> agsub bioturbation.
        if with_tillage
            litter_tillage!(soil, crop)
            tillage_hydraulics!(soil, crop; lpjmlparams = global_params)
        end
        litter_bioturbation!(soil; lpjmlparams = global_params)

        # Preserve LPJmL's albedo/PET-before-snow ordering: radiation uses the
        # snow state present at the start of the day.
        albedo!(pftparameters, crop, soil, pet; maize = maize)  # compute albedo
        petpar!(pet, day_of_year, managed_land.latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables

        snow!(soil, dailyWeather; snowparams = snow_params, lpjmlparams = global_params)

        if water_balance !== nothing
            record_water_balance_after_snow!(water_balance, diagnostic_day, dailyWeather.prec)
        end

        # Thermal properties require current pore volume before phase partitioning.
        pedotransfer!(soil; lpjmlparams = global_params)
        update_surface_litter_properties!(soil; thermalparams = thermal_params)
        soil_temperature!(soil, dailyWeather.temp, climbuf.atemp_mean;
                          thermalparams = thermal_params, snowparams = snow_params)

        # Make same-day mineralization available to crop uptake, following
        # LPJmL's pre-crop litter/SOM and nitrification stage.
        soil_cn_decomposition!(
            soil;
            lpjmlparams = global_params,
            soil_decomp_params = decomp_params,
        )

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)

        annual_output_row = day_of_year == 365 ?
            output_rows.first_annual_row + annual_output_offset : nothing
        harvest_crop!(
            crop, soil, output, managed_land.residue_fraction, day_of_year;
            output_row = output_row,
            annual_output_row = annual_output_row,
        ) # crop harvesting
        route_harvest_residues!(soil, crop)
        annual_output_offset += day_of_year == 365

        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, crop, soil,
                managed_land.residue_fraction,
            )
        end

        # Interception and infiltration precede plant water stress, as in LPJmL.
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec;
                      lpjmlparams = global_params)
        pedotransfer!(soil; lpjmlparams = global_params)
        soil_infiltration!(
            soil,
            crop,
            dailyWeather.prec;
            irrigation = irrigation,
            snowmelt = soil_snow_fluxes(soil).melt,
            air_temperature = dailyWeather.temp,
            lpjmlparams = global_params,
            thermalparams = thermal_params,
        )
        if thermal_balance !== nothing
            record_thermal_balance!(thermal_balance, diagnostic_day, soil)
        end

        if maize
            apar_crop_maize!(pftparameters, crop, pet, soil_snow_prognostic(soil).height) # crop absorbed photosynthetic radiation
        else
            apar_crop!(pftparameters, crop, pet, soil_snow_prognostic(soil).height) # crop absorbed photosynthetic radiation
        end

        temp_stress(pftparameters, pet, crop, dailyWeather.temp;
                    photoparams = photo_params) # temperature stress function

        # C4 photosynthesis
        photosynthesis_C4!(pftparameters, crop, crop_canopy_auxiliary(crop).apar, pet.daylength, dailyWeather.temp;
                           comp_vcmax = true, lpjmlparams = global_params, photoparams = photo_params)

        # Potential conductance at LAMBDA_OPT, followed by the LPJmL
        # water-limited C4 lambda solve on the active CPU/GPU backend.
        transpiration!(crop_fluxes(crop).carbon.water_limited_assimilation, pftparameters, crop, pet, soil, current_co2;
                       lpjmlparams = global_params)
        solve_lambda_c4!(pftparameters, crop, pet, dailyWeather.temp, current_co2;
                         lpjmlparams = global_params, photoparams = photo_params)

        if nitrogen_limit_vcmax
            crop_nitrogen!(crop, pftparameters, soil, crop_photosynthesis_auxiliary(crop).potential_vcmax, dailyWeather.temp;
                           auto_fertilizer = automatic_fertilizer, lpjmlparams = global_params)
            limit_vcmax_by_nitrogen!(crop, pftparameters, dailyWeather.temp;
                                    lpjmlparams = global_params)
        end
        photosynthesis_C4!(pftparameters, crop, crop_canopy_auxiliary(crop).apar, pet.daylength, dailyWeather.temp;
                           comp_vcmax = false, lpjmlparams = global_params, photoparams = photo_params)

        # crop respiration and carbon allocation
        crop_carbon!(
            crop, output, pftparameters, dailyWeather.temp,
            soil_thermal_prognostic(soil).temperature;
            output_row = output_row,
            lpjmlparams = global_params,
        )

        # crop nitrogen allocation
        if nitrogen_limit_vcmax
            allocate_crop_nitrogen!(crop, pftparameters)
        else
            crop_nitrogen!(crop, pftparameters, soil, crop_photosynthesis_auxiliary(crop).vcmax, dailyWeather.temp;
                           auto_fertilizer = automatic_fertilizer,
                           lpjmlparams = global_params) # nitrogen cycle
        end

        evaporation!(pet.eeq, crop, soil; lpjmlparams = global_params)

        # Remove daily plant uptake and soil evaporation after demand/supply calculation.
        soil_evapotranspiration!(soil, crop; irrigation = irrigation)

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
end
