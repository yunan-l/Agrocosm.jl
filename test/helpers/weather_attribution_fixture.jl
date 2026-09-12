function weather_attribution_fixture(cft_id; T = Float64, window_days = 8, phu = 607, sowing_day = 100,
    diurnal_config = nothing, diurnal_amplitude = T(10),
    organ_temperature = false, specific_humidity = T(0.006),
    surface_pressure = T(101325),
    reproductive_sink = false, sterility_rate = nothing,
    sterility_temperature = nothing,
    heat_exposure_config = nothing, daily_statistic_exposure = false,
    terminal_heat = false, filling_rate = nothing, filling_temperature = nothing,
    anthesis_heat = false, cold_sterility = false,
    heat_day_temperature = nothing, heat_day_rate = nothing,
    cold_night_temperature = nothing, cold_night_rate = nothing,
    model_parameters = nothing,
    # Soil texture, so a caller can ask what the column does on something other
    # than the loam this fixture was built on. Defaults are the original values,
    # so every existing caller is unaffected.
    sand = T(0.4), clay = T(0.2), w_sat = T(0.45),
)
    forcing_days = 365 * cld(sowing_day + 200, 365)
    cft = Agrocosm.convert_precision(T, cft_id == 1 ? Agrocosm.cft1 : Agrocosm.cft3)
    # The synthetic fixture runs at a constant 19/25 C, far below any real
    # sterility threshold, so the sink can only be exercised by lowering it here.
    # Terminal heat needs the same treatment and for the same reason: the
    # fixture is far below any real grain-filling threshold too.
    overrides = Dict{Symbol, Any}()
    sterility_rate === nothing || (overrides[:sterility_rate] = T(sterility_rate))
    sterility_temperature === nothing ||
        (overrides[:sterility_temperature] = T(sterility_temperature))
    filling_rate === nothing || (overrides[:filling_rate] = T(filling_rate))
    filling_temperature === nothing ||
        (overrides[:filling_temperature] = T(filling_temperature))
    # The fixture runs at a constant 19/25 C, so the two absolute-threshold
    # mechanisms need their thresholds moved into that range for the same
    # reason the sink and terminal heat do.
    heat_day_temperature === nothing ||
        (overrides[:heat_day_temperature] = T(heat_day_temperature))
    heat_day_rate === nothing || (overrides[:heat_day_rate] = T(heat_day_rate))
    cold_night_temperature === nothing ||
        (overrides[:cold_night_temperature] = T(cold_night_temperature))
    cold_night_rate === nothing ||
        (overrides[:cold_night_rate] = T(cold_night_rate))
    if !isempty(overrides)
        cft = Agrocosm.CFTParameters{T, Int32}(;
            (f => get(overrides, f, getfield(cft, f))
             for f in fieldnames(Agrocosm.CFTParameters))...)
    end
    initial_data = (
        latitude = T[45],
        soilparams = (
            ph = T[6.5], w_sat = fill(T(w_sat), 5, 1),
            sand = reshape(T[sand], 1, 1), clay = reshape(T[clay], 1, 1),
            tdiff_0 = T[0.7], tdiff_15 = T[0.75],
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        ModelState = (
            # Avoid landing exactly on the integer-day harvest threshold in
            # the smooth-gradient fixture. Test that boundary separately.
            crop = (sdate = Int32[sowing_day], phu = T[phu], manure = T[0],
                fertilizer = T[24.55], residuefrac = T[0.67]),
            u0 = (
                swc = reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], 5, 1),
                litc = reshape(T[0.13, 187.5, 225.36], 3, 1),
                fastc = reshape(T[548.97, 368.27, 313.79, 377.55, 344.65], 5, 1),
                slowc = reshape(T[1218.62, 753.33, 660.10, 792.63, 736.38], 5, 1),
                litn = reshape(T[0.0047, 6.47, 9.47], 3, 1),
                fastn = reshape(T[36.60, 24.55, 20.92, 25.17, 22.98], 5, 1),
                slown = reshape(T[81.24, 50.22, 44.01, 52.84, 49.09], 5, 1),
            ),
        ),
    )
    climbuf, crop, pet, soil, management, weather, output =
        Agrocosm.init_states!(cft, initial_data, 1, identity; T)
    climbuf.atemp .= T(10)
    climbuf.temp .= T(10)
    climbuf.atemp_mean .= T(10)
    climate_fields = (
        temp = fill(T(cft_id == 1 ? 19 : 25), forcing_days, 1),
        prec = fill(T(2), forcing_days, 1), sw = fill(T(210), forcing_days, 1),
        lw = fill(T(-40), forcing_days, 1), wind = fill(T(2), forcing_days, 1),
        no3_deposition = fill(T(0.01), forcing_days, 1),
        nh4_deposition = fill(T(0.02), forcing_days, 1), co2 = fill(T(400), cld(forcing_days, 365)),
    )
    # A synthetic, constant diurnal range: this fixture never touches real
    # tasmax/tasmin files, so any nonzero amplitude exercises the sub-daily
    # kernel without needing server data. All THREE exposure sources
    # reconstruct the sub-daily course from it, so the field is added whenever
    # any of them is on; with none on `climate` is exactly as before -- no
    # `diurnal_range` field at all, so existing callers are unaffected.
    needs_range = diurnal_config !== nothing || heat_exposure_config !== nothing ||
        daily_statistic_exposure || anthesis_heat || cold_sterility
    climate = needs_range ? merge(climate_fields, (;
        diurnal_range = fill(T(diurnal_amplitude), forcing_days, 1),
    )) : climate_fields
    # Organ temperature additionally needs humidity and pressure. Constant
    # synthetic values, as with the diurnal range: the point is to exercise the
    # energy balance under AD, not to reproduce a real cell.
    organ_temperature && (climate = merge(climate, (;
        specific_humidity = fill(T(specific_humidity), forcing_days, 1),
        surface_pressure = fill(T(surface_pressure), forcing_days, 1),
    )))
    # A caller-supplied bundle is how a test perturbs one global coefficient,
    # e.g. `senescent_leaf_release`, without editing the defaults.
    parameters = model_parameters === nothing ? Agrocosm.ModelParameters(T) :
        Agrocosm.convert_precision(T, model_parameters)
    state = Agrocosm.model_state(climbuf, crop, pet, soil, management, weather, output)
    processes = Agrocosm.ProcessModules(cft, parameters)
    driver = cft_id == 1 ? Agrocosm.daily_crop_C3! : Agrocosm.daily_crop_C4!
    ordinary = deepcopy(state)
    driver(sowing_day, sowing_day + 200, processes, climate, ordinary;
        fertilizer = :yes, manure = true, with_tillage = true,
        nitrogen_limit_vcmax = true, crop_resp_fix = true, diurnal_config,
        organ_temperature, reproductive_sink,
        heat_exposure_config, daily_statistic_exposure, terminal_heat,
        anthesis_heat, cold_sterility,
        update_vernalization_requirement = false, reuse_output = true)
    events = findall(!iszero, vec(ordinary.output.calendar.harvest_event))
    isempty(events) && error("fixture did not harvest")
    harvest_day = first(events) + sowing_day - 1
    first_day = window_days === :season ? sowing_day + 1 : harvest_day - window_days
    first_day > sowing_day || error("fixture failed before its diagnostic window")
    driver(sowing_day, first_day - 1, processes, climate, state;
        fertilizer = :yes, manure = true, with_tillage = true,
        nitrogen_limit_vcmax = true, crop_resp_fix = true, diurnal_config,
        organ_temperature, reproductive_sink,
        heat_exposure_config, daily_statistic_exposure, terminal_heat,
        anthesis_heat, cold_sterility,
        update_vernalization_requirement = false, reuse_output = true)
    return (; state, cft, parameters, climate, days = first_day:(harvest_day - 1),
        harvest_day, forcing = cat(climate.temp, climate.prec, climate.sw, climate.lw,
            climate.wind; dims = 3))
end
