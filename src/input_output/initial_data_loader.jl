"""
InitialDataLoader(data, data_index, device;
                  load_mineral_nitrogen_restart=false)

Build initial model-state inputs from forcing/parameter datasets. Mineral-N
restart pools are omitted by default because `init_states!` reconstructs NO₃
and NH₄ from slow organic N using the LPJmL fresh-soil initialization rule.
Set `load_mineral_nitrogen_restart=true` only when explicitly restoring a
nitrogen-limited restart state.
"""
function InitialDataLoader(data::NamedTuple,
                           data_index::Vector{Int},
                           device;
                           load_mineral_nitrogen_restart::Bool = false
)


    @unpack latitude, crop, soilparam, initialLPJmL = data

    latitude_set = latitude[data_index] |> device

    crop = (
        sdate = Int32.(crop.sdate[data_index]),
        phu = crop.phu[data_index],
        manure = crop.manure[data_index],
        fertilizer = crop.fertilizer[data_index],
        residuefrac = crop.residuefrac[data_index]
    ) |> device

    soilparam_set = (
        ph = soilparam.soilph[data_index],
        w_sat = soilparam.w_sat[:, data_index],
        sand = reshape(soilparam.sand[data_index], (1, :)),
        clay = reshape(soilparam.clay[data_index], (1, :)),
        # silt = soilparam.silt[data_index],
        tdiff_0 = soilparam.tdiff_0[data_index],
        tdiff_15 = soilparam.tdiff_15[data_index],
        soildepth = soilparam.soildepth,
    ) |> device

    u0_set = (
        swc = initialLPJmL.u0.swc[:, data_index],
        litc = initialLPJmL.u0.litc[:, data_index],
        fastc = initialLPJmL.u0.fastc[:, data_index],
        slowc = initialLPJmL.u0.slowc[:, data_index],
        litn = initialLPJmL.u0.litn[:, data_index],
        fastn = initialLPJmL.u0.fastn[:, data_index],
        slown = initialLPJmL.u0.slown[:, data_index],
    )
    if load_mineral_nitrogen_restart
        u0_set = merge(u0_set, (
            soil_NH4 = initialLPJmL.u0.soil_NH4[:, data_index],
            soil_NO3 = initialLPJmL.u0.soil_NO3[:, data_index],
        ))
    end
    u0_set = u0_set |> device

    model_state = (
        crop = crop,
        c_shift_fast = initialLPJmL.c_shift_fast[:, data_index],
        c_shift_slow = initialLPJmL.c_shift_slow[:, data_index],
        u0 = u0_set
    ) |> device

    InitialData = (
        latitude = latitude_set,
        soilparams = soilparam_set,
        ModelState = model_state
    )

    return InitialData
end
