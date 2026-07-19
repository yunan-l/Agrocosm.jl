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
             mineral_nitrogen_initialization=:restart)

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
                       mineral_nitrogen_initialization::Symbol = :restart
)

    @unpack residue_frac, fastfrac, atmfrac, k_soil10 = lpjmlparams
    @unpack k_litter10, beta_root = PFT

    @unpack latitude, soilparams, ModelState = InitialData

    phu = ModelState.crop.phu
    sdate = ModelState.crop.sdate
    manure = ModelState.crop.manure
    fertilizer = ModelState.crop.fertilizer
    residuefrac = ModelState.crop.residuefrac
    c_shift_fast = ModelState.c_shift_fast
    c_shift_slow = ModelState.c_shift_slow
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
    soil.carbon.shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.carbon.shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.carbon.litter_response = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])

    soil.nitrogen.shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.nitrogen.shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.nitrogen.litter_response = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])

    output = init_output(cell_size, device)

    return climbuf, crop, pet, soil, managed_land, dailyWeather, output
end
