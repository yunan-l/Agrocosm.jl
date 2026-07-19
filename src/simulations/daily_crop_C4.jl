# Purely process-based modelling
"""
daily_crop_C4!(...)

Execute daily forward simulation for C4 crop configuration.
"""
function daily_crop_C4!(day_start, day_end,
                        pftparameters,
                        climate, climbuf, crop, crop_cal, photos, pet, soil, managed_land, 
                        dailyWeather, output;
                        maize = true,
                        irrigation = false,
                        water_balance = nothing
)

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    for day = day_start : day_end

        diagnostic_day = day - day_start + 1

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
        cultivate!(crop, crop_cal, managed_land, soil, day_of_year)

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day) # update climate buffer
        albedo!(pftparameters, crop, pet)  # compute albedo
        petpar!(pet, day_of_year, managed_land.latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, output, managed_land.residuefrac, day_of_year) # crop harvesting
        
        if maize
            apar_crop_maize!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(pftparameters, pet, photos, dailyWeather.temp) # temperature stress function

        # C4 photosynthesis
        photosynthesis_C4!(pftparameters, photos, crop.apar, pet.daylength, dailyWeather.temp; comp_vmax = true)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, output, pftparameters, dailyWeather.temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, pftparameters, soil, photos.vmax, pet.daylength, dailyWeather.temp) # nitrogen cycle         
        
        # evapotranspiration
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        transpiration!(photos.adtmm, pftparameters, crop, pet, soil, dailyWeather.annual_co2)
        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil)

        # soil water cycle
        soil_water!(soil, crop, dailyWeather.prec; irrigation = irrigation)

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, soil, crop)
        end

    end
end
