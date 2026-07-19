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
                        water_balance = nothing,
                        nitrogen_balance = nothing
)

    crop_cal = crop.calendar
    photos = crop.photosynthesis

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    for day = day_start : day_end

        diagnostic_day = day - day_start + 1

        day_of_year = day % 365 != 0 ? day % 365 : 365

        readclimate!(climate, dailyWeather, day)

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

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day) # update climate buffer
        albedo!(pftparameters, crop, pet)  # compute albedo
        petpar!(pet, day_of_year, managed_land.latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)

        harvest_crop!(crop_cal, crop, soil, output, managed_land.residue_fraction, day_of_year) # crop harvesting

        # Interception and infiltration precede plant water stress, as in LPJmL.
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        soil_infiltration!(soil, crop, dailyWeather.prec; irrigation = irrigation)

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
        photosynthesis_C4!(pftparameters, photos, crop.canopy.apar, pet.daylength, dailyWeather.temp; comp_vmax = false)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, output, pftparameters, dailyWeather.temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, pftparameters, soil, photos.vmax, dailyWeather.temp;
                       auto_fertilizer = auto_fertilizer) # nitrogen cycle

        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil; air_temperature = dailyWeather.temp)

        # Remove daily plant uptake and soil evaporation after demand/supply calculation.
        soil_evapotranspiration!(soil, crop; irrigation = irrigation)

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, soil, crop)
        end

        if nitrogen_balance !== nothing
            record_nitrogen_balance_end!(nitrogen_balance, diagnostic_day, crop, soil)
        end

    end
end
