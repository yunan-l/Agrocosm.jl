_pathway_albedo!(::Val{:C3}, cft, crop, soil, pet, maize) =
    albedo!(cft, crop, soil, pet)
_pathway_albedo!(::Val{:C4}, cft, crop, soil, pet, maize) =
    albedo!(cft, crop, soil, pet; maize)

_pathway_apar!(::Val{:C3}, cft, crop, pet, snow_height, maize) =
    apar_crop!(cft, crop, pet, snow_height)
_pathway_apar!(::Val{:C4}, cft, crop, pet, snow_height, maize) =
    apar_crop!(cft, crop, pet, snow_height; maize)

"""Execute one daily crop pathway over the canonical lifecycle state."""
function _daily_crop!(
    pathway::Union{Val{:C3}, Val{:C4}}, start_day, end_day,
    processes::ProcessModules, climate, state::ModelState;
    maize = true,
    irrigation = false,
    manure = false,
    fertilizer = :auto,
    with_tillage = true,
    crop_resp_fix = true,
    nitrogen_limit_vcmax = false,
    diurnal_config = nothing,
    heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
    terminal_heat::Bool = false,
    water_sterility::Bool = false,
    water_filling::Bool = false,
    sowing_mode::Symbol = :prescribed_sdate,
    update_vernalization_requirement::Bool = true,
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
    c_shift_response_sum = nothing,
)
    cftparameters = processes.crop
    model_parameters = processes.global_parameters
    climbuf = state.prognostic.climate
    pet = state.auxiliary.pet
    managed_land = state.inputs.management
    dailyWeather = state.inputs.weather
    output = state.output

    T = eltype(crop_prognostic(state).canopy.lai)
    cftparameters = convert_precision(T, cftparameters)
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

    # Sub-daily integration needs the diurnal temperature range alongside the
    # other forcings. It is read as a row view of the (day, cell) matrix, so no
    # extra buffer and no change to the climate kernels is required.
    # The standalone exposure pass reconstructs the same sub-daily temperature
    # course, so it has the same forcing requirement.
    subdaily_config = diurnal_config === nothing ? heat_exposure_config : diurnal_config
    # The closed form needs the same daily range, but no sub-step machinery.
    exposure_source_count = count((diurnal_config !== nothing,
                                   heat_exposure_config !== nothing,
                                   daily_statistic_exposure))
    if subdaily_config !== nothing || daily_statistic_exposure
        hasproperty(climate, :diurnal_range) || throw(ArgumentError(
            "sub-daily integration requires a `diurnal_range` climate field (tasmax - tasmin)",
        ))
        size(climate.diurnal_range) == size(climate.temp) || throw(DimensionMismatch(
            "diurnal_range must have the same shape as the temperature forcing",
        ))
    end
    # One writer for `heat_exposure_hours`; `SimulationConfiguration` rejects the
    # pair, and this is the guard for callers that bypass it.
    exposure_source_count <= 1 || throw(ArgumentError(
        "sub-daily photosynthesis, the standalone exposure pass and the " *
        "daily-statistic closed form all write heat_exposure_hours; enable at " *
        "most one",
    ))

    # Organ temperature is solved per sub-step, so it has nowhere to live
    # without the sub-daily loop; that combination is rejected here rather than
    # silently ignored. Humidity and pressure ride on `climate` as row views of
    # (day, cell) matrices, exactly like `diurnal_range`.
    # The sterility accumulator is filled per sub-step, so the sink needs the
    # sub-daily loop. It does NOT need organ temperature: without it the
    # accumulator integrates duration at sub-daily air temperature instead of
    # leaf temperature, which is the ablation cell that separates the mechanism
    # from the departure that triggers it. See `runtime_contracts.jl`.
    terminal_heat && exposure_source_count == 0 && throw(ArgumentError(
        "terminal heat requires one of sub-daily photosynthesis, the " *
        "standalone heat-exposure pass or the daily-statistic closed form to " *
        "fill filling_exposure_hours",
    ))
    reproductive_sink && exposure_source_count == 0 && throw(ArgumentError(
        "reproductive sink requires one of sub-daily photosynthesis, the " *
        "standalone heat-exposure pass or the daily-statistic closed form to " *
        "fill heat_exposure_hours",
    ))
    # Unlike the three above, this one needs no exposure field: it reads the
    # forcing directly. What it does need is the daily RANGE, because a daily
    # maximum cannot be recovered from a mean, and a plain daily run is handed a
    # climate with no such field at all.
    anthesis_heat && !hasproperty(climate, :diurnal_range) && throw(ArgumentError(
        "anthesis heat requires a `diurnal_range` climate field (tasmax - tasmin); " *
        "a daily run without it cannot form a daily maximum",
    ))
    # Its mirror, and it needs the range for the mirror reason: a daily MINIMUM
    # cannot be recovered from a mean either.
    cold_sterility && !hasproperty(climate, :diurnal_range) && throw(ArgumentError(
        "cold sterility requires a `diurnal_range` climate field (tasmax - tasmin); " *
        "a daily run without it cannot form a daily minimum",
    ))
    if organ_temperature
        subdaily_config === nothing && throw(ArgumentError(
            "organ temperature requires sub-daily photosynthesis or the standalone " *
            "heat-exposure pass to be enabled",
        ))
        for field in (:specific_humidity, :surface_pressure)
            hasproperty(climate, field) || throw(ArgumentError(
                "organ temperature requires a `$field` climate field",
            ))
            size(getproperty(climate, field)) == size(climate.temp) ||
                throw(DimensionMismatch(
                    "$field must have the same shape as the temperature forcing",
                ))
        end
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
        diurnal = diurnal_config === nothing ? nothing : DiurnalForcing(
            diurnal_config, view(climate.diurnal_range, climate_day, :),
        )
        heat_exposure_forcing = heat_exposure_config === nothing ? nothing :
            DiurnalForcing(
                heat_exposure_config, view(climate.diurnal_range, climate_day, :),
            )
        # The albedo, LAI and conductance arrays are updated in place later in
        # this same day, and this holds references rather than copies, so the
        # kernel reads whatever the day has produced by the time it runs. Only
        # the two climate rows depend on `climate_day` and so are rebuilt here.
        organ = organ_temperature ? OrganTemperatureForcing(
            view(climate.specific_humidity, climate_day, :),
            view(climate.surface_pressure, climate_day, :),
            dailyWeather.wind,
            dailyWeather.swr,
            dailyWeather.lwr,
            pet.albedo,
            crop_canopy_auxiliary(state).actual_lai,
            crop_canopy_auxiliary(state).canopy_conductance,
        ) : nothing

        if carbon_balance !== nothing
            record_carbon_balance_start!(carbon_balance, diagnostic_day, state, state)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_start!(nitrogen_balance, diagnostic_day, state, state)
        end
        if water_balance !== nothing
            record_water_balance_start!(
                water_balance, diagnostic_day, state, dailyWeather.prec,
            )
        end

        # --- Discrete establishment event ---------------------------------
        # Today's climate must enter history before sowing decisions. This is
        # intentionally separate from continuous crop/soil process kernels.
        dynamic_winter_type = isnothing(prescribed_winter_type) ?
            crop_phenology_input(state).winter_type : prescribed_winter_type
        update_climbuf!(
            cftparameters, dailyWeather.temp, climbuf, day;
            prec = dailyWeather.prec,
            dynamic_sowing = sowing_mode === :dynamic_sdate,
            winter_type = dynamic_winter_type,
            update_vernalization_requirement,
        )
        if sowing_mode === :dynamic_sdate
            dynamic_sowing_date!(
                state, climbuf, cftparameters, day_of_year;
                irrigated = irrigation,
                prescribed_winter_type,
            )
        end
        cultivate!(
            state, managed_land, state, day_of_year;
            manure,
            apply_prescribed_fertilizer = fertilizer === :yes,
            defer_second_fertilizer = nitrogen_limit_vcmax,
            prescribed_phu,
            prescribed_winter_type,
            cftparameters = cftparameters,
            lpjmlparams = global_params,
            laimax = cftparameters.laimax,
        )
        if carbon_balance !== nothing
            record_carbon_balance_after_cultivate!(carbon_balance, diagnostic_day, state)
        end

        if with_tillage
            litter_tillage!(state, state)
            tillage_hydraulics!(state, state; lpjmlparams = global_params)
        end
        litter_bioturbation!(state; lpjmlparams = global_params)

        # Radiation uses the snow state present at the start of the day.
        _pathway_albedo!(pathway, cftparameters, state, state, pet, maize)
        petpar!(
            pet, day_of_year, managed_land.latitude, dailyWeather.temp,
            dailyWeather.lwr, dailyWeather.swr,
        )
        sowing_mode === :dynamic_sdate &&
            record_potential_evaporation!(climbuf, pet.eeq, day, global_params)
        # LPJmL keeps melt separate while liquid rain passes through canopy
        # interception; melt joins the throughfall before litter/soil uptake.
        snow!(state, dailyWeather; snowparams = snow_params, lpjmlparams = global_params)

        pedotransfer!(state; lpjmlparams = global_params)
        update_surface_litter_properties!(state; thermalparams = thermal_params)
        soil_temperature!(
            state, dailyWeather.temp, climbuf.atemp_mean;
            thermalparams = thermal_params, snowparams = snow_params,
        )

        # Existing litter/SOM decomposes before today's crop uptake.
        soil_cn_decomposition!(
            state;
            lpjmlparams = global_params,
            soil_decomp_params = decomp_params,
        )
        # Decomposition changes today's litter cover/capacity in both N modes.
        update_surface_litter_properties!(state; thermalparams = thermal_params)
        # Match LPJmL: add atmospheric mineral N after litter/SOM turnover and
        # before the crop's daily nitrogen uptake and nitrogen losses.
        if hasproperty(climate, :no3_deposition) || hasproperty(climate, :nh4_deposition)
            nitrogen_deposition!(
                state, dailyWeather.no3_deposition, dailyWeather.nh4_deposition,
            )
        end
        if nitrogen_limit_vcmax
            # LPJmL applies the second fertilizer/manure dose after today's
            # litter/SOM turnover and atmospheric deposition, before uptake.
            fertilizer!(
                state, managed_land, state, day_of_year;
                fertilizer = fertilizer === :yes,
                manure,
                apply_sowing_dose = false,
                surface_second_manure = nitrogen_limit_vcmax,
                reset_inputs = false,
                lpjmlparams = global_params,
            )
        end
        c_shift_response_sum === nothing ||
            accumulate_c_shift_response!(c_shift_response_sum, state)

        if nitrogen_limit_vcmax
            # LPJmL 5.10 computes `gp_sum` before advancing crop phenology. Keep
            # that raw conductance in today's auxiliary slot; the post-phenology
            # APAR/photosynthesis pass below still supplies the lambda solve.
            _pathway_apar!(
                pathway, cftparameters, state, pet,
                soil_snow_prognostic(state).height, maize,
            )
            temp_stress(
                cftparameters, pet, state, dailyWeather.temp;
                photoparams = photo_params,
            )
            photosynthesis!(
                pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
                pet.daylength, dailyWeather.temp, current_co2, diurnal;
                comp_vcmax = true,
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
            prepare_prephenology_canopy_conductance!(
                cftparameters, state, pet.daylength, current_co2;
                lpjmlparams = global_params,
            )
        end

        # --- Discrete calendar harvest event ------------------------------
        phenology_crop!(
            state, climbuf.V_req, cftparameters, dailyWeather.temp, pet.daylength,
        )
        annual_output_row = day_of_year == 365 ?
            output_rows.first_annual_row + annual_output_offset : nothing
        harvest_crop!(
            state, state, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(state, state)
        annual_output_offset += day_of_year == 365
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, state, state,
                managed_land.residue_fraction,
            )
        end

        interception!(
            state, cftparameters, pet.eeq, dailyWeather.prec;
            lpjmlparams = global_params,
        )
        add_snowmelt_to_precipitation!(
            dailyWeather.prec, soil_snow_fluxes(state).melt,
        )
        if water_balance !== nothing
            record_water_balance_after_snow!(
                water_balance, diagnostic_day, dailyWeather.prec,
            )
        end
        pedotransfer!(state; lpjmlparams = global_params)
        soil_infiltration!(
            state, state, dailyWeather.prec;
            snowmelt = soil_snow_fluxes(state).melt,
            air_temperature = dailyWeather.temp,
            lpjmlparams = global_params,
            thermalparams = thermal_params,
        )
        if thermal_balance !== nothing
            record_thermal_balance!(thermal_balance, diagnostic_day, state)
        end

        _pathway_apar!(
            pathway, cftparameters, state, pet,
            soil_snow_prognostic(state).height, maize,
        )
        temp_stress(
            cftparameters, pet, state, dailyWeather.temp;
            photoparams = photo_params,
        )
        photosynthesis!(
            pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
            pet.daylength, dailyWeather.temp, current_co2, diurnal;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )

        transpiration!(
            crop_fluxes(state).carbon.water_limited_assimilation,
            cftparameters, state, pet, state, current_co2;
            lpjmlparams = global_params,
            use_precomputed_conductance = nitrogen_limit_vcmax,
        )
        solve_lambda!(
            pathway, cftparameters, state, pet, dailyWeather.temp, current_co2;
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        photosynthesis!(
            pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
            pet.daylength, dailyWeather.temp, current_co2, diurnal;
            comp_vcmax = false, organ,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )

        if nitrogen_limit_vcmax
            refresh_potential_canopy_conductance!(
                cftparameters, state, pet.daylength, current_co2,
            )
            acquire_crop_nitrogen!(
                state, cftparameters, state,
                crop_photosynthesis_auxiliary(state).potential_vcmax,
                dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                include_storage_reserve = true,
                biological_fixation = true,
                require_active_photosynthesis = true,
                lpjmlparams = global_params,
            )
            limit_vcmax_by_nitrogen!(
                state, cftparameters, dailyWeather.temp;
                require_active_photosynthesis = true,
                lpjmlparams = global_params,
            )
            photosynthesis!(
                pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
                pet.daylength, dailyWeather.temp, current_co2, diurnal;
                comp_vcmax = false, organ,
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
            recouple_nitrogen_water!(
                pathway, cftparameters, state, pet, state,
                dailyWeather.temp, current_co2;
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
            photosynthesis!(
                pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
                pet.daylength, dailyWeather.temp, current_co2, diurnal;
                comp_vcmax = false, organ,
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
            finalize_nitrogen_limited_transpiration!(
                cftparameters, state, pet, state, current_co2;
                lpjmlparams = global_params,
            )
        end

        # Before allocation: today's exposure must be written, and the harvest
        # index that allocation is about to use has to already reflect any grain
        # set lost. With sub-daily assimilation the exposure arrived with the
        # calls above; with the standalone pass it is taken here, after the
        # canopy state those calls left behind, so both routes integrate the
        # same day's LAI and conductance.
        heat_exposure!(
            cftparameters, state, pet.daylength, dailyWeather.temp,
            heat_exposure_forcing; organ,
        )
        # The range view is only formed when the closed form is on: a daily run
        # is given a climate with no `diurnal_range` field at all, so touching
        # it unconditionally is an error rather than an unused argument.
        if daily_statistic_exposure
            daily_statistic_exposure!(
                cftparameters, state, pet.daylength, dailyWeather.temp,
                view(climate.diurnal_range, climate_day, :), true,
            )
        end
        reproductive_sink && reproductive_sink!(cftparameters, state)
        # Reads the forcing, not an exposure field, so it is independent of which
        # exposure path is on. Same state and same clamp as the sink above, so
        # the two are order-independent.
        anthesis_heat && anthesis_heat!(
            cftparameters, state, dailyWeather.temp,
            view(climate.diurnal_range, climate_day, :),
        )
        # The cold mirror, reading the same two forcing channels and subtracting
        # from the same clamped state, so it is order-independent against all of
        # the above. Placed here rather than beside the vernalization code
        # because this is grain-set damage, not a development requirement.
        cold_sterility && cold_sterility!(
            cftparameters, state, dailyWeather.temp,
            view(climate.diurnal_range, climate_day, :),
        )
        terminal_heat && terminal_heat!(cftparameters, state)
        # Order-independent against the heat sink: both only subtract from
        # `grain_set_fraction` and the state is clamped at zero.
        water_sterility && water_sterility!(cftparameters, state)
        water_filling && water_filling!(cftparameters, state)
        crop_carbon!(
            state, output, cftparameters, dailyWeather.temp,
            soil_thermal_prognostic(state).temperature;
            output_row, crop_resp_fix,
            include_biological_fixation_cost = nitrogen_limit_vcmax,
            lpjmlparams = global_params,
        )
        # --- Discrete failed-crop termination event -----------------------
        # This remains separate from calendar harvest: it is only triggered
        # after today's carbon allocation identifies a non-viable stand.
        terminate_failed_crop!(
            state, state, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(state, state)
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, state, state,
                managed_land.residue_fraction,
            )
        end

        if nitrogen_limit_vcmax
            allocate_crop_nitrogen!(state, cftparameters)
        else
            crop_nitrogen!(
                state, cftparameters, state,
                crop_photosynthesis_auxiliary(state).vcmax, dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                lpjmlparams = global_params,
            )
        end

        evaporation!(
            pet.eeq, state, state;
            lpjmlparams = global_params,
            lpjml_managed_evaporation = true,
        )
        soil_evapotranspiration!(state, state; irrigation)
        record_ecosystem_flux_outputs!(output, state, state; output_row)
        accumulate_season_process_diagnostics!(output, state, state, cftparameters)
        post_crop_nitrogen_losses!(
            state;
            air_temperature = dailyWeather.temp,
            wind_speed = dailyWeather.wind,
            exact_lpjml_volatilization = nitrogen_limit_vcmax,
            lpjmlparams = global_params,
        )

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, state, state)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_end!(
                nitrogen_balance, diagnostic_day, state, state;
                nitrate_deposition = dailyWeather.no3_deposition,
                ammonium_deposition = dailyWeather.nh4_deposition,
            )
        end
        if carbon_balance !== nothing
            record_carbon_balance_end!(carbon_balance, diagnostic_day, state, state)
        end
        # All process kernels for a day share the current backend stream. Wait
        # once at the lifecycle boundary so callers observe a completed daily
        # transition without forcing a host/device barrier after every process.
        synchronize_backend!(crop_prognostic(state).canopy.lai)
    end
    return nothing
end

"""Compatibility entry point for the C3-specialized daily crop transition."""
daily_crop_C3!(args...; kwargs...) = _daily_crop!(Val(:C3), args...; kwargs...)

"""Compatibility entry point for the C4-specialized daily crop transition."""
daily_crop_C4!(args...; kwargs...) = _daily_crop!(Val(:C4), args...; kwargs...)
