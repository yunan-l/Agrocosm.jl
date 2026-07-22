"""
cultivate!(crop, managed_land, soil, day)

Handle sowing-day state transitions and activate crop growth state.
"""
function cultivate_reference!(crop,
                              ml::ManagedLand,
                              soil,
                              day::Int;
                              lpjmlparams::LPJmLParams = lpjmlparams,
                              manure = false,
                              apply_prescribed_fertilizer::Bool = true,
                              laimax = cft1.laimax,
)

    T = eltype(crop_prognostic(crop).canopy.lai)
    sowing = crop_calendar_input(crop).sowing_date .== day
    seed_flaimax = T(0.000083)
    seed_lai = seed_flaimax * T(laimax)
    crop_prognostic(crop).phenology.harvesting .= ifelse.(sowing, false, crop_prognostic(crop).phenology.harvesting)
    crop_events(crop).sowing .= ifelse.(sowing, 1, 0)
    crop_prognostic(crop).phenology.is_growing .= ifelse.(sowing, 1, crop_prognostic(crop).phenology.is_growing)

    # LPJmL allocates a fresh Pftcrop at every cultivation. Agrocosm keeps a
    # persistent struct for GPU execution, so sowing must explicitly reproduce
    # the seasonal part of new_crop() without clearing soil or annual outputs.
    for field in (:vdsum, :husum)
        values = getproperty(crop_prognostic(crop).phenology, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop_phenology_auxiliary(crop).fphu .= ifelse.(
        sowing, zero(T), crop_phenology_auxiliary(crop).fphu,
    )
    crop_prognostic(crop).phenology.senescence .= ifelse.(sowing, false, crop_prognostic(crop).phenology.senescence)
    crop_prognostic(crop).phenology.senescence_previous .= ifelse.(sowing, false, crop_prognostic(crop).phenology.senescence_previous)
    crop_prognostic(crop).phenology.harvesting_previous .= ifelse.(sowing, false, crop_prognostic(crop).phenology.harvesting_previous)
    crop_prognostic(crop).phenology.growing_days .= ifelse.(sowing, 0, crop_prognostic(crop).phenology.growing_days)

    crop_prognostic(crop).canopy.lai .= ifelse.(sowing, seed_lai, crop_prognostic(crop).canopy.lai)
    crop_canopy_auxiliary(crop).flaimax .= ifelse.(sowing, seed_flaimax, crop_canopy_auxiliary(crop).flaimax)
    crop_prognostic(crop).canopy.laimax_adjusted .= ifelse.(sowing, one(T), crop_prognostic(crop).canopy.laimax_adjusted)
    crop_prognostic(crop).canopy.lai_npp_deficit .= ifelse.(sowing, zero(T), crop_prognostic(crop).canopy.lai_npp_deficit)
    crop_prognostic(crop).carbon.biomass .= ifelse.(sowing, T(20), crop_prognostic(crop).carbon.biomass)
    crop_prognostic(crop).carbon.root .= ifelse.(sowing, T(8), crop_prognostic(crop).carbon.root)
    crop_prognostic(crop).carbon.leaf .= ifelse.(sowing, T(0.0113804), crop_prognostic(crop).carbon.leaf)
    crop_prognostic(crop).carbon.storage .= ifelse.(sowing, zero(T), crop_prognostic(crop).carbon.storage)
    crop_prognostic(crop).carbon.pool .= ifelse.(sowing, T(11.9886196), crop_prognostic(crop).carbon.pool)
    init_nitrogen = T(0.7) # C:N ratio of seed = 29
    crop_fluxes(crop).nitrogen.seed_input .= ifelse.(sowing, init_nitrogen, zero(T))
    crop_prognostic(crop).nitrogen.total .= ifelse.(sowing, init_nitrogen, crop_prognostic(crop).nitrogen.total)
    for field in (:leaf, :root, :pool, :storage, :pending_manure,
                  :pending_fertilizer, :stress_sum)
        values = getproperty(crop_prognostic(crop).nitrogen, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop_stress_auxiliary(crop).nitrogen_deficit .= ifelse.(
        sowing, zero(T), crop_stress_auxiliary(crop).nitrogen_deficit,
    )
    crop_prognostic(crop).nitrogen.sufficiency .= ifelse.(sowing, one(T), crop_prognostic(crop).nitrogen.sufficiency)

    for field in (:demand_sum, :supply_sum)
        values = getproperty(crop_prognostic(crop).water, field)
        values .= ifelse.(sowing, zero(T), values)
    end
    crop_stress_auxiliary(crop).water_deficit .= ifelse.(
        sowing, zero(T), crop_stress_auxiliary(crop).water_deficit,
    )
    crop_prognostic(crop).water.sufficiency .= ifelse.(sowing, one(T), crop_prognostic(crop).water.sufficiency)

    # Daily fluxes and diagnostics (NPP, respiration, uptake, demand, current
    # management input, and harvest export) are intentionally not reset here:
    # their owning process overwrites them later in this same daily step.
    fertilizer!(
        crop,
        ml,
        soil,
        day;
        fertilizer = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )

end

function cultivate!(crop,
                    ml::ManagedLand,
                    soil,
                    day::Int;
                    lpjmlparams::LPJmLParams = lpjmlparams,
                    manure = false,
                    apply_prescribed_fertilizer::Bool = true,
                    laimax = cft1.laimax)
    T = eltype(crop_prognostic(crop).canopy.lai)
    launch_1D!(
        cultivate_kernel!,
        crop_calendar_input(crop).sowing_date,
        crop_events(crop).sowing,
        crop_prognostic(crop).phenology.harvesting,
        crop_prognostic(crop).phenology.harvesting_previous,
        crop_prognostic(crop).phenology.is_growing,
        crop_prognostic(crop).phenology.vdsum,
        crop_prognostic(crop).phenology.husum,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.senescence,
        crop_prognostic(crop).phenology.senescence_previous,
        crop_prognostic(crop).phenology.growing_days,
        crop_prognostic(crop).canopy.lai,
        crop_canopy_auxiliary(crop).flaimax,
        crop_prognostic(crop).canopy.laimax_adjusted,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        crop_prognostic(crop).carbon.biomass,
        crop_prognostic(crop).carbon.root,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.pool,
        crop_fluxes(crop).nitrogen.seed_input,
        crop_prognostic(crop).nitrogen.total,
        crop_prognostic(crop).nitrogen.leaf,
        crop_prognostic(crop).nitrogen.root,
        crop_prognostic(crop).nitrogen.pool,
        crop_prognostic(crop).nitrogen.storage,
        crop_prognostic(crop).nitrogen.pending_manure,
        crop_prognostic(crop).nitrogen.pending_fertilizer,
        crop_prognostic(crop).nitrogen.stress_sum,
        crop_prognostic(crop).nitrogen.sufficiency,
        crop_stress_auxiliary(crop).nitrogen_deficit,
        crop_stress_auxiliary(crop).water_deficit,
        crop_prognostic(crop).water.demand_sum,
        crop_prognostic(crop).water.supply_sum,
        crop_prognostic(crop).water.sufficiency,
        T(0.000083),
        T(0.000083) * T(laimax),
        day,
    )
    fertilizer!(
        crop, ml, soil, day;
        fertilizer = apply_prescribed_fertilizer,
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
    else
        seed_input[cell] = zero(T)
    end
end
