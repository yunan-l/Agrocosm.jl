"""
cultivate!(crop, crop_cal, day)

Handle sowing-day state transitions and activate crop growth state.
"""
function cultivate_reference!(crop::Crop,
                              crop_cal::CropCalendar,
                              ml::ManagedLand,
                              soil::Soil,
                              day::Int;
                              lpjmlparams::LPJmLParams = lpjmlparams,
                              manure = false,
                              apply_prescribed_fertilizer::Bool = true,
                              laimax = cft1.laimax,
)

    T = eltype(crop.canopy.lai)
    sowing = crop_cal.sowing_date .== day
    seed_flaimax = T(0.000083)
    seed_lai = seed_flaimax * T(laimax)
    crop.phenology.harvesting .= ifelse.(sowing, false, crop.phenology.harvesting)
    crop_cal.sowing_callback .= ifelse.(sowing, 1, 0)
    crop.phenology.is_growing .= ifelse.(sowing, 1, crop.phenology.is_growing)

    # LPJmL allocates a fresh Pftcrop at every cultivation. Agrocosm keeps a
    # persistent struct for GPU execution, so sowing must explicitly reproduce
    # the seasonal part of new_crop() without clearing soil or annual outputs.
    for field in (:vdsum, :husum, :fphu)
        values = getproperty(crop.phenology, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.phenology.senescence .= ifelse.(sowing, false, crop.phenology.senescence)
    crop.phenology.senescence_previous .= ifelse.(sowing, false, crop.phenology.senescence_previous)
    crop.phenology.harvesting_previous .= ifelse.(sowing, false, crop.phenology.harvesting_previous)
    crop.phenology.growing_days .= ifelse.(sowing, 0, crop.phenology.growing_days)

    crop.canopy.lai .= ifelse.(sowing, seed_lai, crop.canopy.lai)
    crop.canopy.flaimax .= ifelse.(sowing, seed_flaimax, crop.canopy.flaimax)
    crop.canopy.laimax_adjusted .= ifelse.(sowing, one(T), crop.canopy.laimax_adjusted)
    crop.canopy.lai_npp_deficit .= ifelse.(sowing, zero(T), crop.canopy.lai_npp_deficit)
    crop.canopy.phenology_fraction .= ifelse.(sowing, seed_flaimax, crop.canopy.phenology_fraction)
    crop.carbon.biomass .= ifelse.(sowing, T(20), crop.carbon.biomass)
    crop.carbon.root .= ifelse.(sowing, T(8), crop.carbon.root)
    crop.carbon.leaf .= ifelse.(sowing, T(0.0113804), crop.carbon.leaf)
    crop.carbon.storage .= ifelse.(sowing, zero(T), crop.carbon.storage)
    crop.carbon.pool .= ifelse.(sowing, T(11.9886196), crop.carbon.pool)
    init_nitrogen = T(0.7) # C:N ratio of seed = 29
    crop.nitrogen.seed_input .= ifelse.(sowing, init_nitrogen, zero(T))
    crop.nitrogen.total .= ifelse.(sowing, init_nitrogen, crop.nitrogen.total)
    for field in (:leaf, :root, :pool, :storage, :pending_manure,
                  :pending_fertilizer, :stress_sum, :deficit)
        values = getproperty(crop.nitrogen, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.nitrogen.stress .= ifelse.(sowing, one(T), crop.nitrogen.stress)

    for field in (:deficit, :demand_sum, :supply_sum, :waterlogging_days)
        values = getproperty(crop.water, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.water.stress .= ifelse.(sowing, one(T), crop.water.stress)
    crop.water.waterlogging_stress .= ifelse.(sowing, one(T), crop.water.waterlogging_stress)

    # Daily fluxes and diagnostics (NPP, respiration, uptake, demand, current
    # management input, and harvest export) are intentionally not reset here:
    # their owning process overwrites them later in this same daily step.
    fertilizer!(
        crop_cal,
        ml,
        crop,
        soil,
        day;
        enabled = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )

end

function cultivate!(crop::Crop,
                    crop_cal::CropCalendar,
                    ml::ManagedLand,
                    soil::Soil,
                    day::Int;
                    lpjmlparams::LPJmLParams = lpjmlparams,
                    manure = false,
                    apply_prescribed_fertilizer::Bool = true,
                    laimax = cft1.laimax)
    T = eltype(crop.canopy.lai)
    launch_1D!(
        cultivate_kernel!,
        crop_cal.sowing_date,
        crop_cal.sowing_callback,
        crop.phenology.harvesting,
        crop.phenology.harvesting_previous,
        crop.phenology.is_growing,
        crop.phenology.vdsum,
        crop.phenology.husum,
        crop.phenology.fphu,
        crop.phenology.senescence,
        crop.phenology.senescence_previous,
        crop.phenology.growing_days,
        crop.canopy.lai,
        crop.canopy.flaimax,
        crop.canopy.laimax_adjusted,
        crop.canopy.lai_npp_deficit,
        crop.canopy.phenology_fraction,
        crop.carbon.biomass,
        crop.carbon.root,
        crop.carbon.leaf,
        crop.carbon.storage,
        crop.carbon.pool,
        crop.nitrogen.seed_input,
        crop.nitrogen.total,
        crop.nitrogen.leaf,
        crop.nitrogen.root,
        crop.nitrogen.pool,
        crop.nitrogen.storage,
        crop.nitrogen.pending_manure,
        crop.nitrogen.pending_fertilizer,
        crop.nitrogen.stress_sum,
        crop.nitrogen.stress,
        crop.nitrogen.deficit,
        crop.water.deficit,
        crop.water.demand_sum,
        crop.water.supply_sum,
        crop.water.stress,
        crop.water.waterlogging_days,
        crop.water.waterlogging_stress,
        T(0.000083),
        T(0.000083) * T(laimax),
        day,
    )
    fertilizer!(
        crop_cal, ml, crop, soil, day;
        enabled = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function cultivate_kernel!(
    sowing_date::AbstractVector{S},
    sowing_callback::AbstractVector{S},
    harvesting::AbstractVector{B},
    harvesting_previous::AbstractVector{B},
    is_growing::AbstractVector{S},
    vdsum::AbstractVector{T},
    husum::AbstractVector{T},
    fphu::AbstractVector{T},
    senescence::AbstractVector{B},
    senescence_previous::AbstractVector{B},
    growing_days::AbstractVector{S},
    lai::AbstractVector{T},
    flaimax::AbstractVector{T},
    laimax_adjusted::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    phenology_fraction::AbstractVector{T},
    biomass::AbstractVector{T},
    root::AbstractVector{T},
    leaf::AbstractVector{T},
    storage::AbstractVector{T},
    pool::AbstractVector{T},
    seed_input::AbstractVector{T},
    total_nitrogen::AbstractVector{T},
    leaf_nitrogen::AbstractVector{T},
    root_nitrogen::AbstractVector{T},
    pool_nitrogen::AbstractVector{T},
    storage_nitrogen::AbstractVector{T},
    pending_manure::AbstractVector{T},
    pending_fertilizer::AbstractVector{T},
    nitrogen_stress_sum::AbstractVector{T},
    nitrogen_stress::AbstractVector{T},
    nitrogen_deficit::AbstractVector{T},
    water_deficit::AbstractVector{T},
    water_demand_sum::AbstractVector{T},
    water_supply_sum::AbstractVector{T},
    water_stress::AbstractVector{T},
    waterlogging_days::AbstractVector{T},
    waterlogging_stress::AbstractVector{T},
    seed_flaimax::T,
    seed_lai::T,
    day::Integer,
) where {T <: AbstractFloat, S <: Integer, B <: Bool}
    cell = @index(Global)
    sowing = sowing_date[cell] == day
    sowing_callback[cell] = sowing ? one(S) : zero(S)
    if sowing
        harvesting[cell] = false
        harvesting_previous[cell] = false
        is_growing[cell] = one(S)
        vdsum[cell] = zero(T)
        husum[cell] = zero(T)
        fphu[cell] = zero(T)
        senescence[cell] = false
        senescence_previous[cell] = false
        growing_days[cell] = zero(S)
        lai[cell] = seed_lai
        flaimax[cell] = seed_flaimax
        laimax_adjusted[cell] = one(T)
        lai_npp_deficit[cell] = zero(T)
        phenology_fraction[cell] = seed_flaimax
        biomass[cell] = T(20)
        root[cell] = T(8)
        leaf[cell] = T(0.0113804)
        storage[cell] = zero(T)
        pool[cell] = T(11.9886196)
        seed_input[cell] = T(0.7)
        total_nitrogen[cell] = T(0.7)
        leaf_nitrogen[cell] = zero(T)
        root_nitrogen[cell] = zero(T)
        pool_nitrogen[cell] = zero(T)
        storage_nitrogen[cell] = zero(T)
        pending_manure[cell] = zero(T)
        pending_fertilizer[cell] = zero(T)
        nitrogen_stress_sum[cell] = zero(T)
        nitrogen_stress[cell] = one(T)
        nitrogen_deficit[cell] = zero(T)
        water_deficit[cell] = zero(T)
        water_demand_sum[cell] = zero(T)
        water_supply_sum[cell] = zero(T)
        water_stress[cell] = one(T)
        waterlogging_days[cell] = zero(T)
        waterlogging_stress[cell] = one(T)
    else
        seed_input[cell] = zero(T)
    end
end
