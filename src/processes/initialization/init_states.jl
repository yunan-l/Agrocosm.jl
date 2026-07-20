"""
root_distribution(beta_root)

Compute normalized root fractions across the default five soil layers from
the LPJmL-style exponential root profile parameter `beta_root`.
"""
function root_distribution(beta_root::AbstractFloat)

    layerbound = Float32.([200.0, 500.0, 1000.0, 2000.0, 3000.0])

    BOTTOMLAYER = length(layerbound)
    totalroots = 1 - beta_root^(layerbound[BOTTOMLAYER] / 10)
    rootdist = zeros(BOTTOMLAYER)
    rootdist[1] = (1 - beta_root^(layerbound[1] / 10)) / totalroots
    for l in 2:BOTTOMLAYER
        rootdist[l] = (beta_root^(layerbound[l-1] / 10) - beta_root^(layerbound[l] / 10)) / totalroots
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
Crop calendar and photosynthesis state are available as `crop.calendar` and
`crop.photosynthesis`.
"""
function init_states!(PFT::PftParameters,
                       InitialData::NamedTuple,
                       cell_size::Int,
                       device;
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

    dailyWeather = init_weather(cell_size, device)

    climbuf = init_climbuf(cell_size, device)
    crop = init_crop(cell_size, device)
    managed_land = init_managed_land(cell_size, device)
    crop.phenology.phu = copy(phu)
    rootdist = root_distribution(beta_root)
    # idx = crop.phenology.phu .< 0
    # crop.phenology.winter_type[idx] .= true
    # crop.phenology.phu[idx] .= -crop.phenology.phu[idx]
    crop.phenology.winter_type .= ifelse.(crop.phenology.phu .< 0, true, crop.phenology.winter_type)
    crop.phenology.phu .= ifelse.(crop.phenology.phu .< 0, -crop.phenology.phu, crop.phenology.phu)
    crop.water.root_distribution .= device(rootdist)

    crop.calendar.sowing_date = sdate
    managed_land.manure = manure
    managed_land.fertilizer = fertilizer
    managed_land.residue_fraction = residuefrac
    managed_land.latitude = latitude
    pet = init_pet(cell_size, device)
    soil = init_soil(cell_size, soilparams.soildepth, device)
    soil.carbon.litter = copy(u0.litc)
    soil.carbon.fast = copy(u0.fastc)
    soil.carbon.slow = copy(u0.slowc)
    soil.nitrogen.litter = copy(u0.litn)
    soil.nitrogen.fast = copy(u0.fastn)
    soil.nitrogen.slow = copy(u0.slown)
    soil.water.storage = copy(u0.swc)
    initialize_soil_mineral_nitrogen!(
        soil,
        u0,
        mineral_nitrogen_initialization,
    )
    soil.water.saturation_fraction = copy(soilparams.w_sat)
    soil.properties.ph = soilparams.ph
    soil.properties.sand_fraction = soilparams.sand
    soil.properties.clay_fraction = soilparams.clay
    soil.thermal.diffusivity_0 = soilparams.tdiff_0
    soil.thermal.diffusivity_15 = soilparams.tdiff_15

    soil.management.tillage_fraction = device([(1 - residue_frac) 0.0f0 0.0f0; residue_frac 1.0f0 0.0f0; 0.0f0 0.0f0 1.0f0])
    initialize_soil_c_shift!(soil, ModelState, c_shift_initialization)
    days_per_year = 365.0f0
    soil.carbon.litter_response = device(
        [k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    soil.nitrogen.litter_response = device(
        [k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    output = init_output(cell_size, device)

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
