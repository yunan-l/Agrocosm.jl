function weather_attribution_fixture(cft_id; T = Float64, window_days = 8, phu = 607, sowing_day = 100,
    diurnal_config = nothing, diurnal_amplitude = T(10),
)
    forcing_days = 365 * cld(sowing_day + 200, 365)
    cft = Agrocosm.convert_precision(T, cft_id == 1 ? Agrocosm.cft1 : Agrocosm.cft3)
    initial_data = (
        latitude = T[45],
        soilparams = (
            ph = T[6.5], w_sat = fill(T(0.45), 5, 1),
            sand = reshape(T[0.4], 1, 1), clay = reshape(T[0.2], 1, 1),
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
    # kernel without needing server data. `diurnal_config = nothing` (the
    # default) leaves `climate` exactly as before -- no `diurnal_range` field
    # at all, so existing (non-sub-daily) callers of this fixture are unaffected.
    climate = diurnal_config === nothing ? climate_fields : merge(climate_fields, (;
        diurnal_range = fill(T(diurnal_amplitude), forcing_days, 1),
    ))
    parameters = Agrocosm.ModelParameters(T)
    state = Agrocosm.model_state(climbuf, crop, pet, soil, management, weather, output)
    processes = Agrocosm.ProcessModules(cft, parameters)
    driver = cft_id == 1 ? Agrocosm.daily_crop_C3! : Agrocosm.daily_crop_C4!
    ordinary = deepcopy(state)
    driver(sowing_day, sowing_day + 200, processes, climate, ordinary;
        fertilizer = :yes, manure = true, with_tillage = true,
        nitrogen_limit_vcmax = true, crop_resp_fix = true, diurnal_config,
        update_vernalization_requirement = false, reuse_output = true)
    events = findall(!iszero, vec(ordinary.output.calendar.harvest_event))
    isempty(events) && error("fixture did not harvest")
    harvest_day = first(events) + sowing_day - 1
    first_day = window_days === :season ? sowing_day + 1 : harvest_day - window_days
    first_day > sowing_day || error("fixture failed before its diagnostic window")
    driver(sowing_day, first_day - 1, processes, climate, state;
        fertilizer = :yes, manure = true, with_tillage = true,
        nitrogen_limit_vcmax = true, crop_resp_fix = true, diurnal_config,
        update_vernalization_requirement = false, reuse_output = true)
    return (; state, cft, parameters, climate, days = first_day:(harvest_day - 1),
        harvest_day, forcing = cat(climate.temp, climate.prec, climate.sw, climate.lw,
            climate.wind; dims = 3))
end
