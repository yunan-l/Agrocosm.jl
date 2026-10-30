@inline function _enzyme_pathway(cft::Agrocosm.CFTParameters)
    cft.path == 1 && return Val(:C3)
    cft.path == 2 && return Val(:C4)
    throw(ArgumentError("unsupported crop photosynthesis path $(cft.path)"))
end

"""
    _enzyme_continuous_transition!(state, cft, global_parameters, climate, day, observable)

Run the existing daily process kernels between fixed discrete lifecycle events.
Sowing, harvest, failed-crop termination, and output recording are handled by
the outer station driver; this function is the opt-in continuous AD path.
"""
function _enzyme_continuous_transition!(
    state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters,
    global_parameters,
    climate,
    day::Integer,
    observable::Symbol,
    layer_depth,
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
    apply_deferred_prescribed_inputs::Bool = nitrogen_limit_vcmax,
    weather_controls = nothing,
    pathway = Val(:C3),
    lambda_primal_mode = Val(:surrogate),
    irrigation_mode::Symbol = :surface,
)
    irrigation_mode in (:surface, :reset) || throw(ArgumentError(
        "irrigation_mode must be :surface or :reset",
    ))
    T = eltype(Agrocosm.crop_prognostic(state).canopy.lai)
    _enzyme_apply_root_distribution!(state, cft.beta_root)
    model_parameters = _enzyme_model_parameters(T, global_parameters)
    global_params = model_parameters.lpjml
    photo_params = model_parameters.photosynthesis
    snow_params = model_parameters.snow
    thermal_params = model_parameters.soil_thermal
    decomp_params = model_parameters.soil_decomposition
    climbuf = state.prognostic.climate
    pet = state.auxiliary.pet
    managed_land = state.inputs.management
    daily_weather = state.inputs.weather

    current_co2 = Agrocosm.readclimate!(climate, daily_weather, day)
    weather_controls === nothing ||
        Agrocosm.apply_weather_forcing!(daily_weather, weather_controls, day)
    Agrocosm.update_climbuf!(
        cft,
        daily_weather.temp,
        climbuf,
        day;
        update_vernalization_requirement = false,
    )
    # The production driver applies its configured tillage treatment before
    # daily litter mixing. The fixed-event AD path uses the same daily
    # management treatment; sowing and harvest remain outer events.
    Agrocosm.litter_tillage!(state, state)
    Agrocosm.tillage_hydraulics!(state; lpjmlparams = global_params)
    Agrocosm.litter_bioturbation!(state; lpjmlparams = global_params)

    Agrocosm._pathway_albedo!(pathway, cft, state, state, pet, true)
    Agrocosm.petpar!(
        pet,
        day % 365 == 0 ? 365 : day % 365,
        managed_land.latitude,
        daily_weather.temp,
        daily_weather.lwr,
        daily_weather.swr,
    )
    Agrocosm.snow!(state, daily_weather; snowparams = snow_params, lpjmlparams = global_params)
    Agrocosm.pedotransfer!(state; lpjmlparams = global_params)
    Agrocosm.update_surface_litter_properties!(state; thermalparams = thermal_params)
    Agrocosm.soil_temperature!(
        state,
        daily_weather.temp,
        climbuf.atemp_mean;
        thermalparams = thermal_params,
        snowparams = snow_params,
    )
    Agrocosm.soil_cn_decomposition!(
        state;
        lpjmlparams = global_params,
        soil_decomp_params = decomp_params,
    )
    # Keep the post-decomposition litter refresh identical to production.
    Agrocosm.update_surface_litter_properties!(state; thermalparams = thermal_params)
    if hasproperty(climate, :no3_deposition) || hasproperty(climate, :nh4_deposition)
        Agrocosm.nitrogen_deposition!(
            state,
            daily_weather.no3_deposition,
            daily_weather.nh4_deposition,
        )
    end
    if apply_deferred_prescribed_inputs
        # Standard nitrogen-limited trajectories carry the fixed second dose
        # in state; management-adaptation objectives inject their own pulses
        # and explicitly disable this event.
        Agrocosm.fertilizer!(
            state, managed_land, state, weather_controls === nothing ? day : mod1(day, 365);
            fertilizer = true,
            manure = true,
            apply_sowing_dose = false,
            surface_second_manure = true,
            reset_inputs = false,
            lpjmlparams = global_params,
        )
    end
    if nitrogen_limit_vcmax
        # Preserve LPJmL's pre-phenology raw `gp_sum` for the initial water
        # balance. The current-canopy APAR/photosynthesis pass remains below.
        Agrocosm._pathway_apar!(
            pathway, cft,
            state,
            pet,
            Agrocosm.soil_snow_prognostic(state).height,
            true,
        )
        Agrocosm.temp_stress(
            cft,
            pet,
            state,
            daily_weather.temp;
            photoparams = photo_params,
        )
        Agrocosm.photosynthesis!(
            pathway,
            cft,
            state,
            Agrocosm.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            daily_weather.temp,
            current_co2;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        _enzyme_prepare_prephenology_canopy_conductance!(
            state, cft, pet, current_co2, global_params,
        )
    end
    Agrocosm.phenology_crop!(
        state,
        climbuf.V_req,
        cft,
        daily_weather.temp,
        pet.daylength,
    )

    Agrocosm.interception!(
        state,
        cft,
        pet.eeq,
        daily_weather.prec;
        lpjmlparams = global_params,
    )
    Agrocosm.add_snowmelt_to_precipitation!(
        daily_weather.prec,
        Agrocosm.soil_snow_fluxes(state).melt,
    )
    Agrocosm.pedotransfer!(state; lpjmlparams = global_params)
    infiltration_input = irrigation && irrigation_mode === :surface ?
        _enzyme_surface_irrigation_input(
            state, cft, daily_weather.prec,
            Agrocosm.crop_fluxes(state).water.interception,
        ) : daily_weather.prec
    Agrocosm.soil_infiltration!(
        state,
        state,
        infiltration_input;
        snowmelt = Agrocosm.soil_snow_fluxes(state).melt,
        air_temperature = daily_weather.temp,
        lpjmlparams = global_params,
        thermalparams = thermal_params,
    )
    Agrocosm._pathway_apar!(
        pathway, cft,
        state,
        pet,
        Agrocosm.soil_snow_prognostic(state).height,
        true,
    )
    Agrocosm.temp_stress(
        cft,
        pet,
        state,
        daily_weather.temp;
        photoparams = photo_params,
    )
    Agrocosm.photosynthesis!(
        pathway,
        cft,
        state,
        Agrocosm.crop_canopy_auxiliary(state).apar,
        pet.daylength,
        daily_weather.temp,
        current_co2;
        comp_vcmax = true,
        lpjmlparams = global_params,
        photoparams = photo_params,
    )
    Agrocosm.transpiration!(
        Agrocosm.crop_fluxes(state).carbon.water_limited_assimilation,
        cft,
        state,
        pet,
        state,
        current_co2;
        lpjmlparams = global_params,
        use_precomputed_conductance = nitrogen_limit_vcmax,
    )
    _enzyme_solve_lambda_c3!(
        state,
        cft,
        pet,
        daily_weather.temp,
        current_co2,
        global_params,
        photo_params,
        nitrogen_limit_vcmax,
        pathway,
        lambda_primal_mode,
    )
    Agrocosm.photosynthesis!(
        pathway,
        cft,
        state,
        Agrocosm.crop_canopy_auxiliary(state).apar,
        pet.daylength,
        daily_weather.temp,
        current_co2;
        comp_vcmax = false,
        lpjmlparams = global_params,
        photoparams = photo_params,
    )
    if nitrogen_limit_vcmax
        _enzyme_refresh_potential_canopy_conductance!(state, cft, pet, current_co2)
        # Match the nitrogen-limited production path: derive demand and acquire
        # nitrogen from the potential capacity before applying the leaf-N Vcmax
        # constraint to the final photosynthesis evaluation.
        Agrocosm.acquire_crop_nitrogen!(
            state,
            cft,
            state,
            Agrocosm.crop_photosynthesis_auxiliary(state).potential_vcmax,
            daily_weather.temp;
            auto_fertilizer = false,
            include_storage_reserve = true,
            biological_fixation = true,
            require_active_photosynthesis = true,
            lpjmlparams = global_params,
        )
        Agrocosm.limit_vcmax_by_nitrogen!(
            state,
            cft,
            daily_weather.temp;
            require_active_photosynthesis = true,
            lpjmlparams = global_params,
        )
        Agrocosm.photosynthesis!(
            pathway,
            cft,
            state,
            Agrocosm.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            daily_weather.temp,
            current_co2;
            comp_vcmax = false,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        _enzyme_recouple_nitrogen_water_c3!(
            state,
            cft,
            pet,
            daily_weather.temp,
            current_co2,
            global_params,
            photo_params,
            pathway,
            lambda_primal_mode,
        )
        Agrocosm.photosynthesis!(
            pathway,
            cft,
            state,
            Agrocosm.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            daily_weather.temp,
            current_co2;
            comp_vcmax = false,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        _enzyme_finalize_nitrogen_limited_transpiration!(
            state, cft, pet, current_co2, global_params,
        )
    end
    _enzyme_crop_carbon!(
        state,
        cft,
        daily_weather.temp,
        Agrocosm.soil_thermal_prognostic(state).temperature,
        global_params,
        nitrogen_limit_vcmax,
        crop_resp_fix,
    )
    if nitrogen_limit_vcmax
        # Carbon allocation changes organ weights; redistribute the acquired
        # total plant nitrogen exactly where the production driver does.
        Agrocosm.allocate_crop_nitrogen!(state, cft)
    else
        Agrocosm.crop_nitrogen!(
            state,
            cft,
            state,
            Agrocosm.crop_photosynthesis_auxiliary(state).vcmax,
            daily_weather.temp;
            auto_fertilizer = false,
            lpjmlparams = global_params,
        )
    end
    _enzyme_evaporation!(state, pet.eeq, global_params, layer_depth)
    _enzyme_soil_evapotranspiration!(
        state; irrigation = irrigation && irrigation_mode === :reset,
    )
    Agrocosm.post_crop_nitrogen_losses!(
        state;
        air_temperature = daily_weather.temp,
        wind_speed = daily_weather.wind,
        exact_lpjml_volatilization = nitrogen_limit_vcmax,
        lpjmlparams = global_params,
    )
    return nothing
end
