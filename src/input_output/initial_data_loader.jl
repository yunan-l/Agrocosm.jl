"""
InitialDataLoader(data, data_index, device;
                  load_mineral_nitrogen_restart=false,
                  load_c_shift_restart=false)

Build initial model-state inputs from forcing/parameter datasets. Mineral-N
restart pools are omitted by default because `init_states!` reconstructs NO₃
and NH₄ from slow organic N using the LPJmL fresh-soil initialization rule.
Set `load_mineral_nitrogen_restart=true` only when explicitly restoring a
nitrogen-limited restart state.

`c_shift` is also omitted by default. `init_states!` then constructs LPJmL's
fresh-soil distribution internally (0.55 in the top layer and 0.45 shared by
the remaining layers). Set `load_c_shift_restart=true` only when restoring a
post-spin-up or restart distribution.
"""
function InitialDataLoader(data::NamedTuple,
                           data_index::Vector{Int},
                           device;
                           T::Type{<:AbstractFloat} = Float32,
                           load_mineral_nitrogen_restart::Bool = false,
                           load_c_shift_restart::Bool = false
)


    @unpack latitude, crop, soilparam, initialLPJmL = data

    latitude_set = T.(latitude[data_index]) |> device

    crop = (
        sdate = Int32.(crop.sdate[data_index]),
        phu = T.(crop.phu[data_index]),
        manure = T.(crop.manure[data_index]),
        fertilizer = T.(crop.fertilizer[data_index]),
        residuefrac = T.(crop.residuefrac[data_index])
    ) |> device

    soilparam_set = (
        ph = T.(soilparam.soilph[data_index]),
        w_sat = T.(soilparam.w_sat[:, data_index]),
        sand = reshape(T.(soilparam.sand[data_index]), (1, :)),
        clay = reshape(T.(soilparam.clay[data_index]), (1, :)),
        # silt = soilparam.silt[data_index],
        tdiff_0 = T.(soilparam.tdiff_0[data_index]),
        tdiff_15 = T.(soilparam.tdiff_15[data_index]),
        soildepth = T.(soilparam.soildepth),
    ) |> device

    u0_set = (
        swc = T.(initialLPJmL.u0.swc[:, data_index]),
        litc = T.(initialLPJmL.u0.litc[:, data_index]),
        fastc = T.(initialLPJmL.u0.fastc[:, data_index]),
        slowc = T.(initialLPJmL.u0.slowc[:, data_index]),
        litn = T.(initialLPJmL.u0.litn[:, data_index]),
        fastn = T.(initialLPJmL.u0.fastn[:, data_index]),
        slown = T.(initialLPJmL.u0.slown[:, data_index]),
    )
    if load_mineral_nitrogen_restart
        u0_set = merge(u0_set, (
            soil_NH4 = T.(initialLPJmL.u0.soil_NH4[:, data_index]),
            soil_NO3 = T.(initialLPJmL.u0.soil_NO3[:, data_index]),
        ))
    end
    u0_set = u0_set |> device

    model_state = (crop = crop, u0 = u0_set)
    if load_c_shift_restart
        model_state = merge(model_state, (
            c_shift_fast = T.(initialLPJmL.c_shift_fast[:, data_index]),
            c_shift_slow = T.(initialLPJmL.c_shift_slow[:, data_index]),
        ))
    end
    model_state = model_state |> device

    InitialData = (
        latitude = latitude_set,
        soilparams = soilparam_set,
        ModelState = model_state
    )

    return InitialData
end
