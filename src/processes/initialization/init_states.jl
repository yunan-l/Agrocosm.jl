"""
root_distribution(beta_root)

Compute normalized root fractions across the default five soil layers from
the LPJmL-style exponential root profile parameter `beta_root`.
"""
function root_distribution(beta_root::AbstractFloat)
    T = typeof(beta_root)
    layerbound = T[200.0, 500.0, 1000.0, 2000.0, 3000.0]

    BOTTOMLAYER = length(layerbound)
    totalroots = one(T) - beta_root^(layerbound[BOTTOMLAYER] / T(10))
    rootdist = zeros(T, BOTTOMLAYER)
    rootdist[1] = (one(T) - beta_root^(layerbound[1] / T(10))) / totalroots
    for l in 2:BOTTOMLAYER
        rootdist[l] = (
            beta_root^(layerbound[l-1] / T(10)) -
            beta_root^(layerbound[l] / T(10))
        ) / totalroots
    end

    return rootdist
end

"""
initialize_soil_mineral_nitrogen!(soil, u0, strategy)

Initialize soil NO₃ and NH₄ using either restart values or the LPJmL
`initsoil.c` rule. In LPJmL's fresh-soil initialization, each mineral pool in
each layer is initialized to one percent of that layer's slow organic-N pool.
"""
function initialize_soil_mineral_nitrogen!(soil::Soil,
                                           u0::NamedTuple,
                                           strategy::Symbol)
    if strategy === :restart
        if !hasproperty(u0, :soil_NO3) || !hasproperty(u0, :soil_NH4)
            throw(ArgumentError(
                ":restart requires soil_NO3 and soil_NH4; construct inputs " *
                "with load_mineral_nitrogen_restart=true",
            ))
        end
        soil.nitrogen.nitrate .= u0.soil_NO3
        soil.nitrogen.ammonium .= u0.soil_NH4
    elseif strategy === :lpjml_initsoil
        slow_n_fraction = convert(eltype(soil.nitrogen.slow), 0.01)
        soil.nitrogen.nitrate .= soil.nitrogen.slow .* slow_n_fraction
        soil.nitrogen.ammonium .= soil.nitrogen.slow .* slow_n_fraction
    else
        throw(ArgumentError(
            "unknown mineral nitrogen initialization strategy: $strategy; " *
            "use :restart or :lpjml_initsoil",
        ))
    end
    return nothing
end


"""
init_states!(PFT, InitialData, cell_size, device;
             lpjmlparams=lpjmlparams,
             mineral_nitrogen_initialization=:lpjml_initsoil)

Initialize and populate all runtime state structs from static parameters and
input data for one simulation domain.
Returns `(climbuf, crop, pet, soil, managed_land, dailyWeather, output)`.
Crop storage is separated by lifetime into `crop.state`, `crop.fluxes`,
`crop.auxiliary`, and `crop.workspace`.
"""
function init_states!(PFT::PftParameters,
                       InitialData::NamedTuple,
                       cell_size::Int,
                       device;
                       T::Type{<:AbstractFloat} = Float32,
                       lpjmlparams::LPJmLParams = lpjmlparams,
                       mineral_nitrogen_initialization::Symbol = :lpjml_initsoil,
                       c_shift_initialization::Symbol = :lpjml_initsoil
)

    @unpack residue_frac = lpjmlparams
    @unpack k_litter10, beta_root = PFT

    @unpack latitude, soilparams, ModelState = InitialData

    phu = ModelState.crop.phu
    sdate = ModelState.crop.sdate
    manure = ModelState.crop.manure
    fertilizer = ModelState.crop.fertilizer
    residuefrac = ModelState.crop.residuefrac
    u0 = ModelState.u0

    to_float(values) = device(T.(values))
    to_integer(values) = device(Int32.(values))

    dailyWeather = init_weather(T, cell_size, device)
    climbuf = init_climbuf(T, cell_size, device)
    crop = init_crop(T, cell_size, device)
    managed_land = init_managed_land(T, cell_size, device)
    crop.state.phenology.phu = to_float(phu)
    rootdist = root_distribution(T(beta_root))
    # idx = crop.state.phenology.phu .< 0
    # crop.state.phenology.winter_type[idx] .= true
    # crop.state.phenology.phu[idx] .= -crop.state.phenology.phu[idx]
    crop.state.phenology.winter_type .= ifelse.(crop.state.phenology.phu .< 0, true, crop.state.phenology.winter_type)
    crop.state.phenology.phu .= ifelse.(crop.state.phenology.phu .< 0, -crop.state.phenology.phu, crop.state.phenology.phu)
    crop.auxiliary.stress.root_distribution .= device(rootdist)

    crop.state.calendar.sowing_date = to_integer(sdate)
    managed_land.manure = to_float(manure)
    managed_land.fertilizer = to_float(fertilizer)
    managed_land.residue_fraction = to_float(residuefrac)
    managed_land.latitude = to_float(latitude)
    pet = init_pet(T, cell_size, device)
    soil = init_soil(T, cell_size, T.(soilparams.soildepth), device)
    soil.carbon.litter = to_float(u0.litc)
    soil.carbon.fast = to_float(u0.fastc)
    soil.carbon.slow = to_float(u0.slowc)
    soil.nitrogen.litter = to_float(u0.litn)
    soil.nitrogen.fast = to_float(u0.fastn)
    soil.nitrogen.slow = to_float(u0.slown)
    soil.water.storage = to_float(u0.swc)
    initialize_soil_mineral_nitrogen!(
        soil,
        u0,
        mineral_nitrogen_initialization,
    )
    soil.water.saturation_fraction = to_float(soilparams.w_sat)
    soil.properties.ph = to_float(soilparams.ph)
    soil.properties.sand_fraction = to_float(soilparams.sand)
    soil.properties.clay_fraction = to_float(soilparams.clay)
    soil.thermal.diffusivity_0 = to_float(soilparams.tdiff_0)
    soil.thermal.diffusivity_15 = to_float(soilparams.tdiff_15)

    soil.management.tillage_fraction = device(T[
        1 - residue_frac 0 0
        residue_frac 1 0
        0 0 1
    ])
    initialize_soil_c_shift!(soil, ModelState, c_shift_initialization)
    days_per_year = T(365)
    soil.carbon.litter_response = device(
        T[k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    soil.nitrogen.litter_response = device(
        T[k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    output = init_output(T, cell_size, device)

    return climbuf, crop, pet, soil, managed_land, dailyWeather, output
end

"""
    initialize_soil_c_shift!(soil, model_state, strategy)

Initialize the normalized vertical distribution used to route decomposed
litter into fast and slow soil pools.

`:lpjml_initsoil` reproduces LPJmL's fresh-soil initialization: 0.55 in the
top layer and 0.45 distributed uniformly over all remaining layers. `:restart`
restores the fast and slow distributions supplied in `model_state`.
"""
function initialize_soil_c_shift!(soil::Soil,
                                  model_state::NamedTuple,
                                  strategy::Symbol)
    if strategy === :lpjml_initsoil
        layers = size(soil.carbon.shift_fast, 1)
        layers > 1 || throw(ArgumentError("LPJmL c_shift initialization requires at least two soil layers"))
        T = eltype(soil.carbon.shift_fast)
        lower_layer_fraction = T(0.45) / T(layers - 1)

        fill!(soil.carbon.shift_fast, lower_layer_fraction)
        fill!(soil.carbon.shift_slow, lower_layer_fraction)
        soil.carbon.shift_fast[1:1, :] .= T(0.55)
        soil.carbon.shift_slow[1:1, :] .= T(0.55)
    elseif strategy === :restart
        hasproperty(model_state, :c_shift_fast) ||
            throw(ArgumentError("c_shift_initialization=:restart requires ModelState.c_shift_fast"))
        hasproperty(model_state, :c_shift_slow) ||
            throw(ArgumentError("c_shift_initialization=:restart requires ModelState.c_shift_slow"))
        size(model_state.c_shift_fast) == size(soil.carbon.shift_fast) ||
            throw(DimensionMismatch("c_shift_fast must match the soil layer-by-cell shape"))
        size(model_state.c_shift_slow) == size(soil.carbon.shift_slow) ||
            throw(DimensionMismatch("c_shift_slow must match the soil layer-by-cell shape"))
        soil.carbon.shift_fast .= model_state.c_shift_fast
        soil.carbon.shift_slow .= model_state.c_shift_slow
    else
        throw(ArgumentError("unknown c_shift initialization strategy: $strategy"))
    end

    # Carbon and nitrogen use the same fixed routing distribution. Keep it
    # separate from fastfrac and atmfrac, which belong to the daily fluxes.
    soil.nitrogen.shift_fast .= soil.carbon.shift_fast
    soil.nitrogen.shift_slow .= soil.carbon.shift_slow
    return nothing
end
