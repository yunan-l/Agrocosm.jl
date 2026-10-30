function _daily_transition_observable(state::Agrocosm.ModelState, observable::Symbol)
    crop_flux = Agrocosm.crop_fluxes(state)
    soil_carbon_flux = Agrocosm.soil_carbon_fluxes(state)
    soil_water_flux = Agrocosm.soil_water_fluxes(state)
    crop_water_flux = crop_flux.water
    litter_flux = Agrocosm.soil_surface_litter_fluxes(state)
    if observable === :gpp
        return crop_flux.carbon.gross_assimilation[1]
    elseif observable === :reco
        return crop_flux.carbon.respiration[1] +
            crop_flux.carbon.leaf_respiration[1] +
            soil_carbon_flux.heterotrophic_respiration[1]
    elseif observable === :et
        value = crop_water_flux.interception[1] + litter_flux.evaporation[1]
        value += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
        value += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
        value += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
        value += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
        value += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]
        return value
    end
    throw(ArgumentError("unsupported daily transition observable $observable"))
end

"""
    Agrocosm.enzyme_seasonal_loss(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, days, layer_depth, context,
    )

Advance one fixed daily process trajectory and accumulate normalized mean
squared losses for GPP, RECO, and ET. The context masks only decide which
observations contribute; every day still advances the full continuous state.
Discrete events remain outside this objective.
"""
function Agrocosm.enzyme_seasonal_loss(
    theta::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    irrigation::Bool = false,
    irrigation_mode::Symbol = :surface,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    gpp_loss = zero(T)
    reco_loss = zero(T)
    et_loss = zero(T)
    day_index = 1
    while day_index <= length(days)
        _enzyme_continuous_transition!(
            state,
            cft,
            global_parameters,
            climate,
            days[day_index],
            :et,
            layer_depth,
            irrigation,
            nitrogen_limit_vcmax,
            crop_resp_fix,
            true,
            nothing,
            pathway,
            Val(:surrogate),
            irrigation_mode,
        )
        if context.growth_mask[day_index]
            crop_flux = Agrocosm.crop_fluxes(state)
            soil_carbon_flux = Agrocosm.soil_carbon_fluxes(state)
            crop_water_flux = crop_flux.water
            soil_water_flux = Agrocosm.soil_water_fluxes(state)
            litter_flux = Agrocosm.soil_surface_litter_fluxes(state)
            gpp = crop_flux.carbon.gross_assimilation[1]
            reco = crop_flux.carbon.respiration[1] +
                crop_flux.carbon.leaf_respiration[1] +
                soil_carbon_flux.heterotrophic_respiration[1]
            et = crop_water_flux.interception[1] + litter_flux.evaporation[1]
            et += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
            et += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
            et += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
            et += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
            et += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]

            if context.valid_masks.gpp[day_index]
                residual = (gpp - T(context.observations.gpp[day_index])) /
                    T(context.scales.gpp)
                gpp_loss += residual * residual
            end
            if context.valid_masks.reco[day_index]
                residual = (reco - T(context.observations.reco[day_index])) /
                    T(context.scales.reco)
                reco_loss += residual * residual
            end
            if context.valid_masks.et[day_index]
                residual = (et - T(context.observations.et[day_index])) /
                    T(context.scales.et)
                et_loss += residual * residual
            end
        end
        day_index += 1
    end
    gpp_count = T(max(context.counts.gpp, 1))
    reco_count = T(max(context.counts.reco, 1))
    et_count = T(max(context.counts.et, 1))
    return gpp_loss / gpp_count + reco_loss / reco_count + et_loss / et_count
end

mutable struct _ADSeasonalLossAccumulator{T}
    gpp::T
    reco::T
    et::T
end

@inline function _enzyme_add_seasonal_loss!(
    losses::_ADSeasonalLossAccumulator{T},
    state::Agrocosm.ModelState,
    index::Int,
    context::Agrocosm.ADSeasonContext,
) where {T <: AbstractFloat}
    context.growth_mask[index] || return nothing

    crop_flux = Agrocosm.crop_fluxes(state)
    soil_carbon_flux = Agrocosm.soil_carbon_fluxes(state)
    crop_water_flux = crop_flux.water
    soil_water_flux = Agrocosm.soil_water_fluxes(state)
    litter_flux = Agrocosm.soil_surface_litter_fluxes(state)
    gpp = crop_flux.carbon.gross_assimilation[1]
    reco = crop_flux.carbon.respiration[1] +
        crop_flux.carbon.leaf_respiration[1] +
        soil_carbon_flux.heterotrophic_respiration[1]
    et = crop_water_flux.interception[1] + litter_flux.evaporation[1]
    et += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
    et += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
    et += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
    et += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
    et += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]

    if context.valid_masks.gpp[index]
        residual = (gpp - T(context.observations.gpp[index])) / T(context.scales.gpp)
        losses.gpp += residual * residual
    end
    if context.valid_masks.reco[index]
        residual = (reco - T(context.observations.reco[index])) / T(context.scales.reco)
        losses.reco += residual * residual
    end
    if context.valid_masks.et[index]
        residual = (et - T(context.observations.et[index])) / T(context.scales.et)
        losses.et += residual * residual
    end
    return nothing
end

@inline function _enzyme_seasonal_loss_value(
    losses::_ADSeasonalLossAccumulator{T},
    context::Agrocosm.ADSeasonContext,
) where {T <: AbstractFloat}
    gpp_count = T(max(context.counts.gpp, 1))
    reco_count = T(max(context.counts.reco, 1))
    et_count = T(max(context.counts.et, 1))
    return losses.gpp / gpp_count + losses.reco / reco_count + losses.et / et_count
end

"""Run and differentiate one fixed seasonal block.

The context counts remain those of the full season, so summing block losses
reproduces the normalization used by `enzyme_seasonal_loss`.
"""
function _enzyme_seasonal_loss_block(
    theta::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
    crop_resp_fix::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
    irrigation_mode::Symbol = :surface,
) where {T <: AbstractFloat}
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            cft,
            global_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            nitrogen_limit_vcmax,
            crop_resp_fix,
            true,
            nothing,
            pathway,
            Val(:surrogate),
            irrigation_mode,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end
