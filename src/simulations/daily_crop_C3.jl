# Purely process-based modelling
"""
daily_crop_C3!(...)

Execute daily forward simulation for C3 crop configuration.
"""
function daily_crop_C3!(start_day, end_day,
                        pftparameters,
                        climate, climbuf, crop, crop_cal, photos, pet, soil, managed_land, 
                        dailyWeather, output;
                        irrigation = false,
                        auto_fertilizer = true,
                        water_balance = nothing
)

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    for day = start_day : end_day

        diagnostic_day = day - start_day + 1

        day_of_year = day % 365 != 0 ? day % 365 : 365
        
        readclimate!(climate, dailyWeather, day)

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
            apply_prescribed_fertilizer = !auto_fertilizer,
        )

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day) # update climate buffer
        albedo!(pftparameters, crop, pet)  # compute albedo
        petpar!(pet, day_of_year, managed_land.latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, output, managed_land.residuefrac, day_of_year) # crop harvesting

        # Interception and infiltration precede plant water stress, as in LPJmL.
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        soil_infiltration!(soil, crop, dailyWeather.prec; irrigation = irrigation)
        
        apar_crop!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        temp_stress(pftparameters, pet, photos, dailyWeather.temp) # temperature stress function

        # C3 photosynthesis
        photosynthesis_C3!(pftparameters, photos, crop.apar, pet.daylength, dailyWeather.temp, dailyWeather.annual_co2; comp_vmax = true)

        # LPJmL first uses lambda_opt photosynthesis to obtain potential
        # conductance, then constrains conductance by water supply.
        transpiration!(photos.adtmm, pftparameters, crop, pet, soil, dailyWeather.annual_co2)

        # Solve the water-limited lambda on the active backend (CPU or GPU),
        # then recompute photosynthesis with fixed vmax and actual lambda.
        solve_lambda_c3!(pftparameters, photos, crop, pet, dailyWeather.temp, dailyWeather.annual_co2)
        photosynthesis_C3!(pftparameters, photos, crop.apar, pet.daylength, dailyWeather.temp, dailyWeather.annual_co2; comp_vmax = false)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, output, pftparameters, dailyWeather.temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, pftparameters, soil, photos.vmax, dailyWeather.temp;
                       auto_fertilizer = auto_fertilizer) # nitrogen cycle

        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil)

        # Remove daily plant uptake and soil evaporation after demand/supply calculation.
        soil_evapotranspiration!(soil, crop; irrigation = irrigation)

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, soil, crop)
        end

    end
end
