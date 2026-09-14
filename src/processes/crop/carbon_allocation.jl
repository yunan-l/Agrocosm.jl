"""
carbon_allocation!(CFT, crop, photos)

Partition crop biomass among leaf/root/storage/pool carbon compartments.
"""
function carbon_allocation!(CFT::CFTParameters,
                            crop;
                            include_biological_fixation_cost::Bool = false,
                            lpjmlparams::LPJmLParams = lpjmlparams,
)
    # 1D cell-wise allocation; crop_prognostic(crop).carbon.storage provides launch length and kernel arg #1.
    T = eltype(crop_prognostic(crop).carbon.storage)
    kernel_params = (
        FROOTMAX = T(0.4),
        FROOTMIN = T(0.3),
        include_biological_fixation_cost = include_biological_fixation_cost,
        senescent_leaf_release = T(lpjmlparams.senescent_leaf_release),
    )

    launch_1D!(carbon_allocation_kernel!,
               crop_prognostic(crop).carbon.storage,
               crop_events(crop).harvest,
               crop_prognostic(crop).phenology.is_growing,
               crop_prognostic(crop).phenology.growing_days,
               crop_prognostic(crop).nitrogen.stress_sum,
               crop_prognostic(crop).nitrogen.sufficiency,
               crop_stress_auxiliary(crop).nitrogen_deficit,
               crop_stress_auxiliary(crop).water_deficit,
               crop_prognostic(crop).water.sufficiency,
               crop_prognostic(crop).phenology.grain_set_fraction,
               crop_prognostic(crop).phenology.grain_fill_fraction,
               crop_prognostic(crop).phenology.window_assimilate,
               crop_prognostic(crop).phenology.anthesis_reserve,
               crop_prognostic(crop).phenology.filling_progress,
               crop_prognostic(crop).phenology.filling_progress_counted,
               crop_phenology_auxiliary(crop).fphu,
               crop_prognostic(crop).phenology.senescence,
               crop_prognostic(crop).carbon.biomass,
               crop_fluxes(crop).carbon.respiration,
               crop_fluxes(crop).carbon.biological_fixation_cost,
               crop_fluxes(crop).carbon.gross_assimilation,
               crop_fluxes(crop).carbon.leaf_respiration,
               crop_fluxes(crop).carbon.npp,
               crop_prognostic(crop).canopy.lai,
               crop_canopy_auxiliary(crop).actual_lai,
               crop_stress_auxiliary(crop).harvest_index_binding,
               crop_prognostic(crop).carbon.leaf,
               crop_prognostic(crop).carbon.root,
               crop_prognostic(crop).carbon.pool,
               crop_prognostic(crop).canopy.lai_npp_deficit,
               CFT,
               kernel_params)

end

"""Compute daily NPP after leaf, maintenance, and growth respiration."""
@inline compute_crop_npp(
    gross_assimilation::T, leaf_respiration::T, crop_respiration::T,
) where {T <: AbstractFloat} =
    gross_assimilation - leaf_respiration - crop_respiration

"""Compute the seasonal nitrogen sufficiency percentage used by root allocation."""
@inline function compute_seasonal_nitrogen_sufficiency(
    accumulated_sufficiency::T, growing_days::S,
) where {T <: AbstractFloat, S <: Integer}
    return growing_days > zero(S) ? accumulated_sufficiency / T(growing_days) * T(100) : T(100)
end

"""Compute SWAT-style root-carbon fraction from water/N stress and phenology."""
@inline function compute_root_carbon_fraction(
    water_sufficiency::T,
    nitrogen_sufficiency::T,
    phenology_fraction::T,
    root_maximum::T,
    root_minimum::T,
) where {T <: AbstractFloat}
    stress = min(water_sufficiency, nitrogen_sufficiency)
    return root_maximum - (root_minimum * phenology_fraction) * stress /
           (stress + exp(T(6.13) - T(0.0883) * stress))
end

"""Compute water-limited LPJmL harvest index for one crop stand."""
@inline function compute_harvest_index(
    phenology_fraction::T,
    optimal_index::T,
    minimum_index::T,
    water_sufficiency::T,
) where {T <: AbstractFloat}
    potential = T(100) * phenology_fraction /
                (T(100) * phenology_fraction +
                 exp(T(11.1) - T(10) * phenology_fraction))
    optimal = optimal_index > one(T) ? potential * (optimal_index - one(T)) + one(T) :
              potential * optimal_index
    minimum = minimum_index > one(T) ? potential * (minimum_index - one(T)) + one(T) :
              potential * minimum_index
    water_sufficiency >= zero(T) || return optimal
    return (optimal - minimum) * water_sufficiency /
           (water_sufficiency + exp(T(6.13) - T(0.0883) * water_sufficiency)) + minimum
end

"""
    saturating_grain_number(window_assimilate, ceiling, half_carbon)

Grains set per square metre, approaching a genetic ceiling as assimilate during
the critical window stops limiting.

    ceiling * window_assimilate / (window_assimilate + half_carbon)

Two simpler forms were built and measured first, and each failed on real cells.
Linear and unbounded set 8098 maize kernels m-2 against a 3000 target and yielded
19.7 t/ha dry matter - roughly twice the world record - because the coefficient
was derived from the GLOBAL MEAN window NPP while a productive cell overshoots in
proportion. Linear with a hard cap then pinned rice, maize and irrigated wheat at
exactly `ceiling * maximum_grain_carbon` every season, which is a constant yield
and therefore the very defect this mechanism exists to remove.

This form is bounded by construction, has no kink for the reverse pass, and keeps
varying at the top of its range.
"""
@inline function saturating_grain_number(
    window_assimilate::T, ceiling::T, half_carbon::T,
) where {T <: AbstractFloat}
    half_carbon > zero(T) || return zero(T)
    supply = max(window_assimilate, zero(T))
    return ceiling * supply / (supply + half_carbon)
end

"""
    grain_sink_carbon(grain_number, maximum_grain_carbon, fphu, window_end, fill_fraction)

The CERES/DSSAT/APSIM sink: grains already set, times what each can still hold.

Replaces a PRESCRIBED harvest index with an emergent one. `docs/34` measured that
the inherited index is a constant - it reads only season water sufficiency,
through a logistic 97.8% saturated at the lowest value real cells reach, so it
moves through 2% of its own span while binding 83% of maize days. Yield was a
fixed fraction of biomass, which makes every yield loss that is not a biomass
loss - sterility, lodging, sprouting, harvest loss - structurally unrepresentable.

Filling progress runs on `fphu`, which is thermal time, so a hot season completes
filling in fewer DAYS and collects less assimilate: the temperature effect on
grain weight arrives through the source limit rather than through a second
coefficient.

`fill_fraction` is `grain_fill_fraction`, and it multiplies the WEIGHT here while
`grain_set_fraction` multiplies the NUMBER at the call site. That is where this
project's 2x2 was always meant to act; until now both multiplied the constant
index instead, so both were scaling a quantity that carried no information.
"""
@inline function grain_sink_carbon(
    grain_number::T, maximum_grain_carbon::T, progress::T, fill_fraction::T,
) where {T <: AbstractFloat}
    return grain_number * maximum_grain_carbon * progress * fill_fraction
end

"""
    thermal_filling_progress(fphu, window_end)

The unweighted share of grain filling completed, LPJmL's own expression.
"""
@inline function thermal_filling_progress(fphu::T, window_end::T) where {T <: AbstractFloat}
    remaining = one(T) - window_end
    return remaining > zero(T) ?
           clamp((fphu - window_end) / remaining, zero(T), one(T)) : one(T)
end

"""
    filling_weight(water_sufficiency, exponent)

How much of a day's filling a grain actually deposits, given that day's water
sufficiency. `exponent = 0` returns one and makes the weighted integral the
thermal progress bitwise.
"""
@inline function filling_weight(water_sufficiency::T, exponent::T) where {T <: AbstractFloat}
    exponent > zero(T) || return one(T)
    return clamp(water_sufficiency, zero(T), one(T))^exponent
end

"""Compute and mass-cap storage carbon after leaf/root allocation."""
@inline function compute_storage_carbon(
    biomass::T,
    leaf_carbon::T,
    root_carbon::T,
    root_fraction::T,
    harvest_index::T,
    optimal_index::T,
) where {T <: AbstractFloat}
    leaf_carbon + root_carbon < biomass || return zero(T)
    candidate = optimal_index > one(T) ?
                (one(T) - one(T) / harvest_index) * (one(T) - root_fraction) * biomass :
                harvest_index * (one(T) - root_fraction) * biomass
    return min(candidate, biomass - leaf_carbon - root_carbon)
end

"""True when the harvest index, and not the carbon mass cap or the grain already
deposited, set storage carbon today.

The expressions are duplicated from `compute_storage_carbon` verbatim rather than
refactored out of it, so that function stays bitwise untouched - it is on the
differentiated path and every ablation rung's equivalence rests on it. The branch
selector is `optimal_index`, NOT `harvest_index`: writing `harvest_index > one(T)`
would flip the formula on every cell where a reproductive mechanism pushed the
index below one, which is exactly the population this diagnostic exists to count.
"""
@inline function harvest_index_binds(
    biomass::T,
    leaf_carbon::T,
    root_carbon::T,
    root_fraction::T,
    harvest_index::T,
    optimal_index::T,
    deposited::T,
) where {T <: AbstractFloat}
    leaf_carbon + root_carbon < biomass || return false
    candidate = optimal_index > one(T) ?
                (one(T) - one(T) / harvest_index) * (one(T) - root_fraction) * biomass :
                harvest_index * (one(T) - root_fraction) * biomass
    mass_cap = biomass - leaf_carbon - root_carbon
    return candidate < mass_cap && candidate > deposited
end

@kernel inbounds = true function carbon_allocation_kernel!(
                                           crop_stoc::AbstractArray{T},
                                           crop_harvest::AbstractArray{S},
                                           crop_isgrowing::AbstractArray{S},
                                           crop_growingdays::AbstractArray{S},
                                           crop_vscal_sum::AbstractArray{T},
                                           crop_vscal::AbstractArray{T},
                                           crop_ndf::AbstractArray{T},
                                           crop_wdf::AbstractArray{T},
                                           crop_wscal::AbstractArray{T},
                                           crop_grain_set::AbstractArray{T},
                                           crop_grain_fill::AbstractArray{T},
                                           crop_window_assimilate::AbstractArray{T},
                                           crop_anthesis_reserve::AbstractArray{T},
                                           crop_filling_progress::AbstractArray{T},
                                           crop_filling_counted::AbstractArray{T},
                                           crop_fphu::AbstractArray{T},
                                           crop_senescence::AbstractArray{B},
                                           crop_biomass::AbstractArray{T},
                                           crop_resp::AbstractArray{T},
                                           crop_bnf_cost::AbstractArray{T},
                                           photos_agd::AbstractArray{T},
                                           photos_rd::AbstractArray{T},
                                           crop_npp::AbstractArray{T},
                                           crop_lai::AbstractArray{T},
                                           crop_actual_lai::AbstractArray{T},
                                           crop_hi_binding::AbstractArray{T},
                                           crop_leafc::AbstractArray{T},
                                           crop_rootc::AbstractArray{T},
                                           crop_poolc::AbstractArray{T},
                                           crop_lai_nppdeficit::AbstractArray{T},
                                           CFT::CFTParameters,
                                           kernel_params
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)

    @unpack sla, hiopt, himin = CFT
    @unpack grain_number_half_carbon, maximum_grain_carbon, maximum_grain_number = CFT
    @unpack reserve_remobilisation, filling_stress_exponent = CFT
    @unpack flowering_start, flowering_end = CFT
    @unpack FROOTMAX, FROOTMIN, include_biological_fixation_cost = kernel_params
    @unpack senescent_leaf_release = kernel_params

    # Diagnostic only, and daily-owned: zeroed every day so a flag cannot survive
    # from yesterday on a cell that took no allocation path today.
    crop_hi_binding[cell] = zero(T)

    if crop_isgrowing[cell] == 1
        # LPJmL preserves the potential phenological LAI and applies the NPP
        # deficit only when actual LAI is consumed or reported.
        actual_lai = max(zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell])
        # Complete crop carbon cost: leaf respiration plus maintenance/growth
        # respiration, including root respiration.
        crop_npp[cell] = compute_crop_npp(
            photos_agd[cell], photos_rd[cell], crop_resp[cell],
        )
        if include_biological_fixation_cost
            crop_npp[cell] -= crop_bnf_cost[cell]
        end
        if ((crop_biomass[cell] + crop_npp[cell]) <= T(0.0001)) || ((actual_lai <= zero(T)) && (!crop_senescence[cell]))
            # LPJmL reports `negbm` here. The daily driver then harvests the
            # remaining pools and removes the failed crop stand.
            crop_poolc[cell] += crop_npp[cell]
            crop_biomass[cell] += crop_npp[cell]
            crop_harvest[cell] = one(S)
        else
            crop_biomass[cell] += crop_npp[cell]
            crop_vscal_sum[cell] += crop_vscal[cell]
            crop_ndf[cell] = compute_seasonal_nitrogen_sufficiency(
                crop_vscal_sum[cell], crop_growingdays[cell],
            )

            # Root carbon follows SWAT-style stress-scaled partitioning.
            froot = compute_root_carbon_fraction(
                crop_wdf[cell], crop_ndf[cell], crop_fphu[cell], FROOTMAX, FROOTMIN,
            )
            crop_rootc[cell] = froot * crop_biomass[cell]

            # Grain filling is irreversible: carbon already deposited in the
            # storage organ cannot be taken back to build leaves. Reserve it
            # before leaf allocation, capped by the above-ground carbon that
            # actually exists.
            #
            # Without this, the two phenological branches disagree about
            # priority. Senescence (below) protects storage and makes leaves
            # give way; pre-senescence used to do the opposite, letting leaves
            # claim all above-ground carbon and driving storage to zero via
            # `compute_storage_carbon`'s `leaf + root < biomass` guard. Measured
            # on the Michigan soybean cell, that reversed 6.7 gC of already-set
            # grain across four days with the canopy still standing. Crops that
            # never approach carbon limitation (wheat, rice, maize here) never
            # took that branch, so this changes nothing for them.
            deposited = min(crop_stoc[cell],
                            max(zero(T), crop_biomass[cell] - crop_rootc[cell]))
            leaf_available = max(zero(T),
                                 crop_biomass[cell] - crop_rootc[cell] - deposited)

            # Leaf carbon is constrained by LAI and SLA; in senescence it is mass-balanced.
            if !crop_senescence[cell]
                if leaf_available >= (crop_lai[cell] / sla)
                    crop_leafc[cell] = crop_lai[cell] / sla
                    crop_lai_nppdeficit[cell] = zero(T)
                else
                    crop_leafc[cell] = leaf_available
                    crop_lai_nppdeficit[cell] = crop_lai[cell] - crop_leafc[cell] * sla
                end
            else
                # Senesced leaf area no longer stands, so the carbon that
                # supported it stops being leaf carbon. LPJmL keeps `leaf`
                # frozen at its last pre-senescence value and only ever trims it
                # through the mass-balance clamp below, which leaves a crop
                # holding a canopy's worth of carbon at LAI = 0 while
                # `compute_storage_carbon`'s `biomass - leaf - root` cap denies
                # that same carbon to the grain. The released carbon goes to the
                # mobile pool via the balance closure further down, where the
                # harvest index still decides how much of it becomes grain, so
                # this removes a bookkeeping constraint rather than adding a
                # flux. `senescent_leaf_release = 0` is LPJmL, bitwise.
                standing_leaf = max(
                    zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell],
                ) / sla
                surplus = crop_leafc[cell] - standing_leaf
                # A branch rather than `leafc -= release * max(0, surplus)`.
                # Once the canopy is gone both sides are zero, and the
                # subtract-a-max form does arithmetic on that exact tie every
                # remaining day; the branch touches nothing there.
                if surplus > zero(T)
                    crop_leafc[cell] -= senescent_leaf_release * surplus
                end
                if (crop_leafc[cell] + crop_rootc[cell] + crop_stoc[cell]) > crop_biomass[cell]
                    crop_leafc[cell] = crop_biomass[cell] - crop_rootc[cell] - crop_stoc[cell]
                end
                if crop_leafc[cell] < zero(T)
                    crop_leafc[cell] = zero(T)
                end
            end

            # Storage carbon (harvest index branch) is computed after leaf/root partitioning.
            # Grain that failed to set caps the harvest index, and grain that
            # filled poorly caps it again. The carbon denied to storage stays in
            # the pool below, so biomass is conserved while yield falls -- the
            # signature of reproductive heat damage that a photosynthesis-only
            # path cannot produce.
            #
            # The two factors MULTIPLY because they limit different things in
            # sequence: how many grains were set, and how far each of those
            # filled. Adding them would let one mechanism repair the other's
            # damage, and taking the minimum would make the less severe of the
            # two free.
            # Grain number accumulates from assimilate supply inside the critical
            # window and is fixed once the window closes - the defining property
            # of the CERES structure, and the reason a poor flowering fortnight
            # caps yield however good the rest of the season is.
            # The state accumulates window ASSIMILATE; the grain number is the
            # saturating function of it, evaluated where it is used. Storing the
            # driver rather than the result keeps the response function in one
            # place and lets it change without a state migration.
            if T(grain_number_half_carbon) > zero(T) &&
               crop_fphu[cell] > T(flowering_start) && crop_fphu[cell] < T(flowering_end)
                crop_window_assimilate[cell] += max(zero(T), crop_npp[cell])
            end
            # Stem carbon standing at anthesis, recorded once. What the grain may
            # take of it is `reserve_remobilisation`; everything the canopy fixes
            # afterwards is available in full.
            above_ground = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell]
            if crop_fphu[cell] >= T(flowering_start) &&
               crop_anthesis_reserve[cell] <= zero(T)
                crop_anthesis_reserve[cell] = max(above_ground, zero(T))
            end
            # At `reserve_remobilisation = 1` this is `above_ground` exactly,
            # which is the inherited cap bitwise.
            fillable = above_ground - (one(T) - T(reserve_remobilisation)) *
                       min(crop_anthesis_reserve[cell], above_ground)
            # Each day's filling weighted by that day's water sufficiency. At
            # exponent 0 the weight is one and this is the thermal progress
            # bitwise, which is the ablation contract.
            thermal_progress = thermal_filling_progress(crop_fphu[cell], T(flowering_end))
            increment = max(thermal_progress - crop_filling_counted[cell], zero(T))
            crop_filling_progress[cell] += increment *
                filling_weight(crop_wscal[cell], T(filling_stress_exponent))
            crop_filling_counted[cell] = thermal_progress
            effective_progress = crop_filling_progress[cell]
            hi = compute_harvest_index(crop_fphu[cell], T(hiopt), T(himin), crop_wdf[cell]) *
                 crop_grain_set[cell] * crop_grain_fill[cell]
            # Never below what is already deposited: the harvest-index formula
            # describes how much grain the crop is filling towards, not a
            # quantity that can be un-filled. Mass still closes, because
            # `deposited` was capped at the available above-ground carbon and
            # leaves were allocated from what remained after it.
            # `grain_number_half_carbon = 0` takes the branch below and leaves the
            # inherited index bitwise untouched, which is the ablation contract
            # every mechanism in this project ships with.
            sink = T(grain_number_half_carbon) > zero(T) ?
                min(grain_sink_carbon(
                        saturating_grain_number(
                            crop_window_assimilate[cell], T(maximum_grain_number),
                            T(grain_number_half_carbon)) * crop_grain_set[cell],
                        T(maximum_grain_carbon), effective_progress,
                        crop_grain_fill[cell]),
                    fillable) :
                compute_storage_carbon(
                    crop_biomass[cell], crop_leafc[cell], crop_rootc[cell], froot, hi,
                    T(hiopt))
            crop_stoc[cell] = max(deposited, sink)
            # Counts the days the harvest index actually bound. A mechanism that
            # multiplies `hi` changes nothing on a day the mass cap or the
            # already-deposited grain bound instead, so this is what decides
            # whether such a mechanism can matter at all.
            crop_hi_binding[cell] = harvest_index_binds(
                crop_biomass[cell], crop_leafc[cell], crop_rootc[cell], froot, hi,
                T(hiopt), deposited,
            ) ? one(T) : zero(T)

            # Pool carbon closes biomass balance and is clipped during senescence if negative.
            crop_poolc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell] - crop_stoc[cell]
            # pool can become negative during senescence
            if crop_senescence[cell] && crop_poolc[cell] < zero(T)
                if (crop_stoc[cell] + crop_poolc[cell]) < zero(T)
                    crop_poolc[cell] += crop_stoc[cell]
                    crop_stoc[cell] = zero(T)
                    if (crop_rootc[cell] + crop_poolc[cell]) < zero(T)
                        crop_poolc[cell] += crop_rootc[cell]
                        crop_rootc[cell] = zero(T) # remainder negative pool must be compensated by leaves,
                        crop_leafc[cell] += crop_poolc[cell]
                        crop_poolc[cell] = zero(T)
                    else
                        crop_rootc[cell] += crop_poolc[cell]
                        crop_poolc[cell] = zero(T)
                    end
                else
                    crop_stoc[cell] += crop_poolc[cell]
                    crop_poolc[cell] = zero(T)
                end
            end
        end

    else
        crop_leafc[cell] = zero(T)
        crop_rootc[cell] = zero(T)
        crop_stoc[cell] = zero(T)
        crop_poolc[cell] = zero(T)
        crop_npp[cell] = zero(T)
        crop_biomass[cell] = zero(T)
        crop_vscal_sum[cell] = zero(T)
        crop_ndf[cell] = zero(T)
        crop_lai_nppdeficit[cell] = zero(T)
    end

    crop_actual_lai[cell] = max(zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell])

end
