"""
transpiration!(photos_adtmm, CFT, crop, pet, soil, co2; lpjmlparams=lpjmlparams)

Compute water demand/supply balance and layer-resolved transpiration uptake.
"""
function transpiration!(photos_adtmm::AbstractArray{T},
                        CFT::CFTParameters,
                        crop,
                        pet::PetPar,
                        soil,
                        co2::AbstractArray{T};
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        use_precomputed_conductance::Bool = false,
) where {T <: AbstractFloat}

    # Root-zone weighted water availability is accumulated inside the cell
    # kernel, avoiding a separate broadcast and reduction array every day.
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop_prognostic(crop).carbon.root))
    # demand = ifelse.(crop_canopy_auxiliary(crop).canopy_conductance .> 0, (1 .- crop_canopy_auxiliary(crop).canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop_canopy_auxiliary(crop).canopy_conductance), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    kernel_params = (
        lpjmlparams = lpjmlparams,
        soil_layers = 5,
        use_precomputed_conductance = use_precomputed_conductance,
    )

    launch_1D!(water_demand_supply_kernel!,
               crop_canopy_auxiliary(crop).canopy_conductance,
               photos_adtmm,
               co2,
               pet.daylength,
               crop_canopy_auxiliary(crop).fpar,
               crop_fluxes(crop).water.transpiration_layer,
               crop_prognostic(crop).water.demand_sum,
               crop_prognostic(crop).water.supply_sum,
               crop_stress_auxiliary(crop).water_deficit,
               crop_prognostic(crop).water.sufficiency,
               crop_prognostic(crop).carbon.root,
               crop_canopy_auxiliary(crop).canopy_wet,
               crop_prognostic(crop).phenology.is_growing,
               pet.eeq,
               crop_root_input(crop).distribution,
               crop_root_auxiliary(crop).zone_available_water,
               soil_water_auxiliary(soil).relative_content,
               soil_water_auxiliary(soil).holding_capacity_storage,
               CFT,
               kernel_params)

end

"""
    prepare_prephenology_canopy_conductance!(CFT, crop, daylength, co2)

Store LPJmL's raw `gp_sum` crop conductance after photosynthesis has been
evaluated with the canopy state present before today's phenology update.
"""
function prepare_prephenology_canopy_conductance!(
    CFT::CFTParameters,
    crop,
    daylength::AbstractArray{T},
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    launch_1D!(
        prepare_prephenology_canopy_conductance_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        co2,
        daylength,
        crop_canopy_auxiliary(crop).fpar,
        T(CFT.gmin),
        T(lpjmlparams.LAMBDA_OPT),
    )
    return nothing
end

@kernel inbounds = true function prepare_prephenology_canopy_conductance_kernel!(
    conductance::AbstractArray{T},
    assimilation::AbstractArray{T},
    co2::AbstractArray{T},
    daylength::AbstractArray{T},
    fpar::AbstractArray{T},
    gmin::T,
    lambda_optimum::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    co2_cell = co2[length(co2) == 1 ? 1 : cell]
    conductance[cell] = compute_canopy_conductance(
        assimilation[cell], co2_cell, daylength[cell], fpar[cell], gmin,
        lambda_optimum,
    )
end

"""Return whether LPJmL repeats the water-limited lambda solve after N stress."""
@inline nitrogen_water_recoupling_required(
    potential_conductance::T,
    nitrogen_limited_conductance::T,
    demand::T,
    supply::T,
) where {T <: AbstractFloat} =
    potential_conductance - nitrogen_limited_conductance > T(0.01) &&
    demand - supply > T(0.1)

"""Store LPJmL's pre-N-limitation `gc_new` from potential assimilation."""
function refresh_potential_canopy_conductance!(
    CFT::CFTParameters,
    crop,
    daylength::AbstractArray{T},
    co2::AbstractArray{T},
) where {T <: AbstractFloat}
    launch_1D!(
        refresh_potential_canopy_conductance_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).fpar,
        crop_prognostic(crop).phenology.is_growing,
        daylength,
        co2,
        T(CFT.gmin),
    )
    return nothing
end

@kernel inbounds = true function refresh_potential_canopy_conductance_kernel!(
    conductance::AbstractArray{T},
    potential_assimilation::AbstractArray{T},
    lambda::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    fpar::AbstractArray{T},
    is_growing::AbstractArray{S},
    daylength::AbstractArray{T},
    co2::AbstractArray{T},
    gmin::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        conductance[cell] = compute_canopy_conductance(
            potential_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, lambda[cell],
        )
    end
end

"""
    recouple_nitrogen_water!(pathway, CFT, crop, pet, soil, temperature, co2)

Apply LPJmL's conditional second lambda solve after leaf nitrogen has limited
Vcmax. The first water pass remains responsible for seasonal water-stress
diagnostics; this correction only updates today's conductance and lambda.
"""
function recouple_nitrogen_water!(
    pathway::Union{Val{:C3}, Val{:C4}},
    CFT::CFTParameters,
    crop,
    pet::PetPar,
    soil,
    temperature::AbstractArray{T},
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
    photoparams::PhotoParams = photoparams,
) where {T <: AbstractFloat}
    kernel_params = (
        pathway = pathway,
        b = T(CFT.b),
        emax = T(CFT.emax),
        fpc = T(CFT.fpc),
        gmin = T(CFT.gmin),
        soil_layers = 5,
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )
    launch_1D!(
        nitrogen_water_recoupling_kernel!,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).apar,
        crop_canopy_auxiliary(crop).canopy_wet,
        crop_prognostic(crop).phenology.is_growing,
        crop_prognostic(crop).carbon.root,
        crop_root_input(crop).distribution,
        soil_water_auxiliary(soil).relative_content,
        pet.daylength,
        pet.eeq,
        temperature,
        co2,
        kernel_params,
    )
    return nothing
end

@kernel inbounds = true function nitrogen_water_recoupling_kernel!(
    lambda::AbstractArray{T},
    vcmax::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    limited_assimilation::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    canopy_wet::AbstractArray{T},
    is_growing::AbstractArray{S},
    root_carbon::AbstractArray{T},
    root_distribution::AbstractArray{T},
    soil_water::AbstractArray{M},
    daylength::AbstractArray{T},
    equilibrium_evaporation::AbstractArray{T},
    temperature::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    @unpack pathway, b, emax, fpc, gmin, soil_layers, lpjmlparams, photoparams = kernel_params
    @unpack ALPHAM, GM = lpjmlparams

    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        previous_lambda = lambda[cell]
        potential_conductance = conductance[cell]
        limited_conductance = compute_canopy_conductance(
            limited_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, previous_lambda,
        )
        demand = compute_transpiration_demand(
            canopy_wet[cell], equilibrium_evaporation[cell], T(ALPHAM), T(GM),
            limited_conductance,
        )
        root_water = zero(T)
        for layer in 1:soil_layers
            root_water += soil_water[layer, cell] * root_distribution[layer]
        end
        supply = compute_transpiration_supply(emax, root_water, root_carbon[cell]) * fpc
        conductance[cell] = limited_conductance

        if nitrogen_water_recoupling_required(
            potential_conductance, limited_conductance, demand, supply,
        )
            constrained_conductance = compute_actual_canopy_conductance(
                limited_conductance, supply, demand, canopy_wet[cell],
                equilibrium_evaporation[cell], T(ALPHAM), T(GM),
            )
            gpd, fac = compute_canopy_water_supply(
                daylength[cell], constrained_conductance, gmin, fpar[cell], co2_cell,
            )
            if gpd > T(1e-5) && daylength[cell] > zero(T) && co2_cell > zero(T)
                if pathway isa Val{:C3}
                    lambda[cell] = compute_lambda_c3_solution(
                        fac, vcmax[cell], temperature_stress[cell], b, co2_cell,
                        temperature[cell], apar[cell], daylength[cell], lpjmlparams,
                        photoparams, previous_lambda, 20,
                    )
                else
                    lambda[cell] = compute_lambda_c4_solution(
                        fac, vcmax[cell], temperature_stress[cell], b,
                        temperature[cell], apar[cell], daylength[cell], lpjmlparams,
                        photoparams, previous_lambda, 20,
                    )
                end
            end
        end
    end
end

"""Overwrite today's transpiration layers using final N-limited photosynthesis."""
function finalize_nitrogen_limited_transpiration!(
    CFT::CFTParameters,
    crop,
    pet::PetPar,
    soil,
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    kernel_params = (gmin = T(CFT.gmin), soil_layers = 5,
                     lpjmlparams = lpjmlparams)
    launch_1D!(
        finalize_nitrogen_limited_transpiration_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).water.transpiration_layer,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).canopy_wet,
        crop_prognostic(crop).phenology.is_growing,
        crop_root_input(crop).distribution,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).holding_capacity_storage,
        pet.daylength,
        pet.eeq,
        co2,
        kernel_params,
    )
    return nothing
end

@kernel inbounds = true function finalize_nitrogen_limited_transpiration_kernel!(
    conductance::AbstractArray{T},
    transpiration_layer::AbstractArray{T},
    limited_assimilation::AbstractArray{T},
    lambda::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    fpar::AbstractArray{T},
    canopy_wet::AbstractArray{T},
    is_growing::AbstractArray{S},
    root_distribution::AbstractArray{T},
    soil_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    daylength::AbstractArray{T},
    equilibrium_evaporation::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    @unpack gmin, soil_layers, lpjmlparams = kernel_params
    @unpack ALPHAM, GM = lpjmlparams

    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        final_conductance = compute_canopy_conductance(
            limited_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, lambda[cell],
        )
        demand = compute_transpiration_demand(
            canopy_wet[cell], equilibrium_evaporation[cell], T(ALPHAM), T(GM),
            final_conductance,
        )
        root_water = zero(T)
        for layer in 1:soil_layers
            root_water += soil_water[layer, cell] * root_distribution[layer]
        end
        transpiration = root_water > zero(T) ? demand * fpar[cell] / root_water : zero(T)
        for layer in 1:soil_layers
            transpiration_layer[layer, cell], _ = compute_layer_transpiration(
                transpiration, root_distribution[layer], soil_water[layer, cell],
                holding_storage[layer, cell],
            )
        end
        conductance[cell] = final_conductance
    end
end

"""Compute LPJmL canopy conductance from water-limited assimilation."""
@inline function compute_canopy_conductance(
    water_limited_assimilation::T,
    co2::T,
    daylength::T,
    fpar::T,
    minimum_conductance::T,
    lambda_optimum::T,
) where {T <: AbstractFloat}
    co2_bar = co2 * T(1e-5)
    co2_bar > zero(T) && daylength > zero(T) || return zero(T)
    denominator = co2_bar * (one(T) - lambda_optimum) * hour2sec(daylength)
    return T(1.6) * water_limited_assimilation / denominator +
           minimum_conductance * fpar
end

"""Compute root-water-limited transpiration supply for one crop column."""
@inline compute_transpiration_supply(
    maximum_supply::T, root_water::T, root_carbon::T,
) where {T <: AbstractFloat} =
    maximum_supply * root_water * (one(T) - exp(T(-0.04) * root_carbon))

"""Compute canopy transpiration demand from conductance and wet-canopy fraction."""
@inline function compute_transpiration_demand(
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
    conductance::T,
) where {T <: AbstractFloat}
    conductance > zero(T) || return zero(T)
    return (one(T) - canopy_wet) * equilibrium_evaporation * alpha /
           (one(T) + (conductance_shape * alpha) / conductance)
end

"""Return the LPJmL 0–100 seasonal water-sufficiency diagnostic."""
@inline function compute_water_sufficiency(
    supplied::T, demanded::T,
) where {T <: AbstractFloat}
    demanded > zero(T) || return T(100)
    return clamp(T(100) * supplied / demanded, zero(T), T(100))
end

"""Cap one layer's transpiration extraction by its plant-available water."""
@inline function compute_layer_transpiration(
    transpiration::T,
    root_fraction::T,
    relative_water::T,
    holding_storage::T,
) where {T <: AbstractFloat}
    unconstrained = transpiration * root_fraction * relative_water
    capacity = relative_water * holding_storage
    return min(unconstrained, capacity), unconstrained > capacity
end

"""Recover actual canopy conductance after layer-wise uptake capping."""
@inline function compute_actual_canopy_conductance(
    current_conductance::T,
    actual_supply::T,
    demand::T,
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
) where {T <: AbstractFloat}
    actual_supply < demand && equilibrium_evaporation > zero(T) || return current_conductance
    denominator = (one(T) - canopy_wet) * equilibrium_evaporation * alpha - actual_supply
    return denominator > zero(T) ? conductance_shape * alpha * actual_supply / denominator : zero(T)
end

@kernel inbounds = true function water_demand_supply_kernel!(
                                             crop_gp::AbstractArray{T},
                                             photos_adtmm::AbstractArray{T},
                                             co2::AbstractArray{T},
                                             daylength::AbstractArray{T},
                                             crop_fpar::AbstractArray{T},
                                             crop_trans_layer::AbstractArray{T},
                                             crop_w_demandsum::AbstractArray{T},
                                             crop_w_supplysum::AbstractArray{T},
                                             crop_wdf::AbstractArray{T},
                                             crop_wscal::AbstractArray{T},
                                             crop_rootc::AbstractArray{T},
                                             crop_canopy_wet::AbstractArray{T},
                                             crop_isgrowing::AbstractArray{S},
                                             pet_eeq::AbstractArray{T},
                                             crop_rootdist::AbstractArray{T},
                                             crop_rootzone_available_water::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             CFT::CFTParameters,
                                             kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers, use_precomputed_conductance = kernel_params
    @unpack ALPHAM, GM, LAMBDA_OPT = lpjmlparams
    @unpack fpc, emax, gmin = CFT

    co2_index = length(co2) == 1 ? 1 : cell
    if !use_precomputed_conductance
        crop_gp[cell] = compute_canopy_conductance(
            photos_adtmm[cell], co2[co2_index], daylength[cell], crop_fpar[cell],
            T(gmin), T(LAMBDA_OPT),
        )
    end

    wr = zero(T)
    rootzone_water = zero(T)
    for l in 1:soil_layers
        wr += soil_w[l, cell] * crop_rootdist[l]
        if l <= 3
            rootzone_water += soil_w[l, cell] * soil_whcs[l, cell] * crop_rootdist[l]
        end
    end
    crop_rootzone_available_water[cell] = rootzone_water

    if crop_isgrowing[cell] == 1
        supply = compute_transpiration_supply(T(emax), wr, crop_rootc[cell])
        demand = compute_transpiration_demand(
            crop_canopy_wet[cell], pet_eeq[cell], T(ALPHAM), T(GM), crop_gp[cell],
        )

        crop_w_demandsum[cell] += demand
        if supply > demand
            crop_w_supplysum[cell] += demand
        else
            crop_w_supplysum[cell] += supply
        end

        crop_wdf[cell] = compute_water_sufficiency(
            crop_w_supplysum[cell], crop_w_demandsum[cell],
        )

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            crop_wscal[cell] = (emax * wr) / (pet_eeq[cell] * ALPHAM / (one(T) + (GM * ALPHAM) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        # Potential transpiration constrained by demand/supply and canopy fraction.
        if wr > 0
            transp = min(supply, demand) / wr * fpc
        else
            transp = zero(T)
        end

        transp_cor = zero(T)

        # Apply layer-wise extraction cap so uptake does not exceed layer storage.
        if transp > 0
            for l in 1:soil_layers
                transp_tmp, capped = compute_layer_transpiration(
                    transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
                )
                transp_cor += transp_tmp
                capped && transp_cor < T(1e-5) && (transp_cor = zero(T))
            end
        else
            transp_cor = zero(T)
        end

        if wr > 0
            transp = transp_cor / wr
        else
            transp = zero(T)
        end

        # LPJmL recomputes actual canopy conductance after layer extraction.
        # Store it in `gp`; downstream lambda solving consumes this actual value.
        actual_supply = fpc > zero(T) ? transp_cor / fpc : zero(T)
        crop_gp[cell] = compute_actual_canopy_conductance(
            crop_gp[cell], actual_supply, demand, crop_canopy_wet[cell], pet_eeq[cell],
            T(ALPHAM), T(GM),
        )

        # Distribute corrected transpiration back to layers by root distribution.
        for l in 1:soil_layers
            crop_trans_layer[l, cell], _ = compute_layer_transpiration(
                transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
            )
        end
    else
        crop_gp[cell] = zero(T)
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = zero(T)
        end
        crop_w_demandsum[cell] = zero(T)
        crop_w_supplysum[cell] = zero(T)
        crop_wdf[cell] = zero(T)
        # Neutral stress for an absent stand; is_growing still gates all fluxes.
        crop_wscal[cell] = one(T)
    end
end
