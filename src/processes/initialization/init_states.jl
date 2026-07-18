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
init_states!(PFT, InitialData, cell_size, device; lpjmlparams=lpjmlparams)

Initialize and populate all runtime state structs from static parameters and
input data for one simulation domain.
Returns `(climbuf, crop, crop_cal, photos, pet, soil, managed_land, dailyWeather, output)`.
"""
function init_states!(PFT::PftParameters,
                       InitialData::NamedTuple,
                       cell_size::Int,
                       device;
                       lpjmlparams::LPJmLParams = lpjmlparams
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
    crop, crop_cal, managed_land, photos = init_crop(cell_size, device)
    crop.phu = copy(phu)
    rootdist = root_distribution(beta_root)
    # idx = crop.phu .< 0
    # crop.wtype[idx] .= true
    # crop.phu[idx] .= -crop.phu[idx]
    crop.wtype .= ifelse.(crop.phu .< 0, true, crop.wtype)
    crop.phu .= ifelse.(crop.phu .< 0, -crop.phu, crop.phu)
    crop.rootdist .= device(rootdist)

    crop_cal.sdate = sdate
    managed_land.manure = manure
    managed_land.fertilizer = fertilizer
    managed_land.residuefrac = residuefrac
    managed_land.latitude = latitude
    pet = init_pet(cell_size, device)
    soil = init_soil(cell_size, soilparams.soildepth, device)
    soil.litc = copy(u0.litc)
    soil.fastc = copy(u0.fastc)
    soil.slowc = copy(u0.slowc)
    soil.litn = copy(u0.litn)
    soil.fastn = copy(u0.fastn)
    soil.slown = copy(u0.slown)
    soil.swc = copy(u0.swc)
    soil.NO3 = copy(u0.soil_NO3)
    soil.NH4 = copy(u0.soil_NH4)
    soil.wsat = copy(soilparams.w_sat)
    soil.ph = soilparams.ph
    soil.sand = soilparams.sand
    soil.clay = soilparams.clay
    soil.tdiff_0 = soilparams.tdiff_0
    soil.tdiff_15 = soilparams.tdiff_15

    soil.tillage_frac = device([(1 - residue_frac) 0.0f0 0.0f0; residue_frac 1.0f0 0.0f0; 0.0f0 0.0f0 1.0f0])
    soil.c_shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.c_shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.response_litc = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])

    soil.n_shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.n_shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.response_litn = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])
    
    output = init_output(cell_size, device)

    return climbuf, crop, crop_cal, photos, pet, soil, managed_land, dailyWeather, output
end
