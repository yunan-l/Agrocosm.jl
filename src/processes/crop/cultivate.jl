"""
cultivate!(crop, managed_land, soil, day)

Handle sowing-day state transitions and activate crop growth state.
"""
function cultivate_reference!(crop::Crop,
                              ml::ManagedLand,
                              soil::Soil,
                              day::Int;
                              lpjmlparams::LPJmLParams = lpjmlparams,
                              manure = false,
                              apply_prescribed_fertilizer::Bool = true,
                              laimax = cft1.laimax,
)

    T = eltype(crop.state.canopy.lai)
    sowing = crop.state.calendar.sowing_date .== day
    seed_flaimax = T(0.000083)
    seed_lai = seed_flaimax * T(laimax)
    crop.state.phenology.harvesting .= ifelse.(sowing, false, crop.state.phenology.harvesting)
    crop.events.sowing .= ifelse.(sowing, 1, 0)
    crop.state.phenology.is_growing .= ifelse.(sowing, 1, crop.state.phenology.is_growing)

    # LPJmL allocates a fresh Pftcrop at every cultivation. Agrocosm keeps a
    # persistent struct for GPU execution, so sowing must explicitly reproduce
    # the seasonal part of new_crop() without clearing soil or annual outputs.
    for field in (:vdsum, :husum, :fphu)
        values = getproperty(crop.state.phenology, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.state.phenology.senescence .= ifelse.(sowing, false, crop.state.phenology.senescence)
    crop.state.phenology.senescence_previous .= ifelse.(sowing, false, crop.state.phenology.senescence_previous)
    crop.state.phenology.harvesting_previous .= ifelse.(sowing, false, crop.state.phenology.harvesting_previous)
    crop.state.phenology.growing_days .= ifelse.(sowing, 0, crop.state.phenology.growing_days)

    crop.state.canopy.lai .= ifelse.(sowing, seed_lai, crop.state.canopy.lai)
    crop.auxiliary.canopy.flaimax .= ifelse.(sowing, seed_flaimax, crop.auxiliary.canopy.flaimax)
    crop.state.canopy.laimax_adjusted .= ifelse.(sowing, one(T), crop.state.canopy.laimax_adjusted)
    crop.state.canopy.lai_npp_deficit .= ifelse.(sowing, zero(T), crop.state.canopy.lai_npp_deficit)
    crop.auxiliary.canopy.phenology_fraction .= ifelse.(sowing, seed_flaimax, crop.auxiliary.canopy.phenology_fraction)
    crop.state.carbon.biomass .= ifelse.(sowing, T(20), crop.state.carbon.biomass)
    crop.state.carbon.root .= ifelse.(sowing, T(8), crop.state.carbon.root)
    crop.state.carbon.leaf .= ifelse.(sowing, T(0.0113804), crop.state.carbon.leaf)
    crop.state.carbon.storage .= ifelse.(sowing, zero(T), crop.state.carbon.storage)
    crop.state.carbon.pool .= ifelse.(sowing, T(11.9886196), crop.state.carbon.pool)
    init_nitrogen = T(0.7) # C:N ratio of seed = 29
    crop.fluxes.nitrogen.seed_input .= ifelse.(sowing, init_nitrogen, zero(T))
    crop.state.nitrogen.total .= ifelse.(sowing, init_nitrogen, crop.state.nitrogen.total)
    for field in (:leaf, :root, :pool, :storage, :pending_manure,
                  :pending_fertilizer, :stress_sum)
        values = getproperty(crop.state.nitrogen, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.auxiliary.stress.nitrogen_deficit .= ifelse.(
        sowing, zero(T), crop.auxiliary.stress.nitrogen_deficit,
    )
    crop.auxiliary.stress.nitrogen .= ifelse.(sowing, one(T), crop.auxiliary.stress.nitrogen)

    for field in (:demand_sum, :supply_sum, :waterlogging_days)
        values = getproperty(crop.state.water, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop.auxiliary.stress.water_deficit .= ifelse.(
        sowing, zero(T), crop.auxiliary.stress.water_deficit,
    )
    crop.auxiliary.stress.water .= ifelse.(sowing, one(T), crop.auxiliary.stress.water)
    crop.auxiliary.stress.waterlogging .= ifelse.(sowing, one(T), crop.auxiliary.stress.waterlogging)

    # Daily fluxes and diagnostics (NPP, respiration, uptake, demand, current
    # management input, and harvest export) are intentionally not reset here:
    # their owning process overwrites them later in this same daily step.
    fertilizer!(
        crop,
        ml,
        soil,
        day;
        enabled = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )

end

function cultivate!(crop::Crop,
                    ml::ManagedLand,
                    soil::Soil,
                    day::Int;
                    lpjmlparams::LPJmLParams = lpjmlparams,
                    manure = false,
                    apply_prescribed_fertilizer::Bool = true,
                    laimax = cft1.laimax)
    T = eltype(crop.state.canopy.lai)
    launch_1D!(
        cultivate_kernel!,
        crop.state.calendar.sowing_date,
        crop.events.sowing,
        crop.state.phenology.harvesting,
        crop.state.phenology.harvesting_previous,
        crop.state.phenology.is_growing,
        crop.state.phenology.vdsum,
        crop.state.phenology.husum,
        crop.state.phenology.fphu,
        crop.state.phenology.senescence,
        crop.state.phenology.senescence_previous,
        crop.state.phenology.growing_days,
        crop.state.canopy.lai,
        crop.auxiliary.canopy.flaimax,
        crop.state.canopy.laimax_adjusted,
        crop.state.canopy.lai_npp_deficit,
        crop.auxiliary.canopy.phenology_fraction,
        crop.state.carbon.biomass,
        crop.state.carbon.root,
        crop.state.carbon.leaf,
        crop.state.carbon.storage,
        crop.state.carbon.pool,
        crop.fluxes.nitrogen.seed_input,
        crop.state.nitrogen.total,
        crop.state.nitrogen.leaf,
        crop.state.nitrogen.root,
        crop.state.nitrogen.pool,
        crop.state.nitrogen.storage,
        crop.state.nitrogen.pending_manure,
        crop.state.nitrogen.pending_fertilizer,
        crop.state.nitrogen.stress_sum,
        crop.auxiliary.stress.nitrogen,
        crop.auxiliary.stress.nitrogen_deficit,
        crop.auxiliary.stress.water_deficit,
        crop.state.water.demand_sum,
        crop.state.water.supply_sum,
        crop.auxiliary.stress.water,
        crop.state.water.waterlogging_days,
        crop.auxiliary.stress.waterlogging,
        T(0.000083),
        T(0.000083) * T(laimax),
        day,
    )
    fertilizer!(
        crop, ml, soil, day;
        enabled = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function cultivate_kernel!(
    sowing_date::AbstractVector{S},
    sowing_event::AbstractVector{S},
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
    sowing_event[cell] = sowing ? one(S) : zero(S)
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
