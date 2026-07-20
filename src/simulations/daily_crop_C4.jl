# Purely process-based modelling
"""
daily_crop_C4!(...)

Execute daily forward simulation for C4 crop configuration.
"""
function daily_crop_C4!(day_start, day_end,
                        pftparameters,
                        climate, climbuf, crop, pet, soil, managed_land,
                        dailyWeather, output;
                        maize = true,
                        irrigation = false,
                        manure = false,
                        auto_fertilizer = true,
                        nitrogen_limit_vmax = false,
                        water_balance = nothing,
                        nitrogen_balance = nothing,
                        carbon_balance = nothing,
                        thermal_balance = nothing
)

    crop_cal = crop.calendar
    photos = crop.photosynthesis

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    annual_rows = count(day -> day % 365 == 0, day_start:day_end)
    output_rows = prepare_output_block!(output, day_end - day_start + 1, annual_rows)
    annual_output_offset = 0

    for day = day_start : day_end

        diagnostic_day = day - day_start + 1
        output_row = output_rows.first_daily_row + diagnostic_day - 1

        day_of_year = day % 365 != 0 ? day % 365 : 365

        readclimate!(climate, dailyWeather, day)

        if carbon_balance !== nothing
            record_carbon_balance_start!(carbon_balance, diagnostic_day, crop, soil)
        end

        if nitrogen_balance !== nothing
            record_nitrogen_balance_start!(nitrogen_balance, diagnostic_day, crop, soil)
        end

        if water_balance !== nothing
            record_water_balance_start!(water_balance, diagnostic_day, soil, dailyWeather.prec)
        end

        # snow
        snow!(soil, dailyWeather)

        if water_balance !== nothing
            record_water_balance_after_snow!(water_balance, diagnostic_day, dailyWeather.prec)
        end

        # initial crop variables in sowing day and fertilizer
        cultivate!(
            crop,
            crop_cal,
            managed_land,
            soil,
            day_of_year;
            manure = manure,
            apply_prescribed_fertilizer = !auto_fertilizer,
        )

        if carbon_balance !== nothing
            record_carbon_balance_after_cultivate!(carbon_balance, diagnostic_day, crop)
        end

        # LPJmL tills existing litter at cultivation, then applies the daily
        # agtop -> agsub bioturbation transfer before surface-litter physics.
        litter_tillage!(soil, crop_cal)
        litter_bioturbation!(soil)

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day) # update climate buffer
        albedo!(pftparameters, crop, pet)  # compute albedo
        petpar!(pet, day_of_year, managed_land.latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        update_surface_litter_properties!(soil)
        # Thermal properties require current pore volume before phase partitioning.
        pedotransfer!(soil)
        soil_temperature!(soil, dailyWeather.temp, climbuf.atemp_mean)

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)

        annual_output_row = day_of_year == 365 ?
            output_rows.first_annual_row + annual_output_offset : nothing
        harvest_crop!(
            crop_cal, crop, soil, output, managed_land.residue_fraction, day_of_year;
            output_row = output_row,
            annual_output_row = annual_output_row,
        ) # crop harvesting
        annual_output_offset += day_of_year == 365

        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, crop, soil,
                managed_land.residue_fraction,
            )
        end

        # Interception and infiltration precede plant water stress, as in LPJmL.
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        soil_infiltration!(
            soil,
            crop,
            dailyWeather.prec;
            irrigation = irrigation,
            snowmelt = soil.snow.melt,
            air_temperature = dailyWeather.temp,
        )
        if thermal_balance !== nothing
            record_thermal_balance!(thermal_balance, diagnostic_day, soil)
        end

        if maize
            apar_crop_maize!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(pftparameters, pet, photos, dailyWeather.temp) # temperature stress function

        # C4 photosynthesis
        photosynthesis_C4!(pftparameters, photos, crop.canopy.apar, pet.daylength, dailyWeather.temp; comp_vmax = true)

        # Potential conductance at LAMBDA_OPT, followed by the LPJmL
        # water-limited C4 lambda solve on the active CPU/GPU backend.
        transpiration!(photos.water_limited_assimilation, pftparameters, crop, pet, soil, dailyWeather.annual_co2)
        solve_lambda_c4!(pftparameters, photos, crop, pet, dailyWeather.temp, dailyWeather.annual_co2)

        if nitrogen_limit_vmax
            crop_nitrogen!(crop, pftparameters, soil, photos.potential_vmax, dailyWeather.temp;
                           auto_fertilizer = auto_fertilizer)
            limit_vmax_by_nitrogen!(crop, pftparameters, dailyWeather.temp)
        end
        photosynthesis_C4!(pftparameters, photos, crop.canopy.apar, pet.daylength, dailyWeather.temp; comp_vmax = false)

        # crop respiration and carbon allocation
        crop_carbon!(
            photos, crop, output, pftparameters, dailyWeather.temp;
            output_row = output_row,
        )

        # crop nitrogen allocation
        if nitrogen_limit_vmax
            allocate_crop_nitrogen!(crop, pftparameters)
        else
            crop_nitrogen!(crop, pftparameters, soil, photos.vmax, dailyWeather.temp;
                           auto_fertilizer = auto_fertilizer) # nitrogen cycle
        end

        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil;
                       air_temperature = dailyWeather.temp,
                       wind_speed = dailyWeather.wind)

        # Remove daily plant uptake and soil evaporation after demand/supply calculation.
        soil_evapotranspiration!(soil, crop; irrigation = irrigation)

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
