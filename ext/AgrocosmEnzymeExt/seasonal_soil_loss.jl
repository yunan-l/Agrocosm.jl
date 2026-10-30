"""
    Agrocosm.enzyme_seasonal_soil_loss(
        theta_soil, state, base_cft, base_parameters, soil_parameter_names,
        climate, days, layer_depth, context,
    )

Accumulate the fixed-event seasonal loss while differentiating only selected
soil fields inside `ModelParameters`. The CFT and all non-soil parameter
bundles remain constant.
"""
function Agrocosm.enzyme_seasonal_soil_loss(
    theta_soil::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    irrigation::Bool = false,
    irrigation_mode::Symbol = :surface,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in eachindex(days)
        _enzyme_continuous_transition!(
            state,
            base_cft,
            model_parameters,
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

"""
    Agrocosm.enzyme_seasonal_joint_loss(
        theta_cft, theta_soil, state, base_cft, base_parameters,
        cft_parameter_names, soil_parameter_names, climate, days,
        layer_depth, context,
    )

Accumulate one fixed-event seasonal loss while differentiating selected CFT
and soil fields together. Existing CFT-only and soil-only entry points remain
unchanged.
"""
function Agrocosm.enzyme_seasonal_joint_loss(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    irrigation::Bool = false,
    irrigation_mode::Symbol = :surface,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    cft = _replace_cft_parameters(base_cft, theta_cft, cft_parameter_names)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in eachindex(days)
        _enzyme_continuous_transition!(
            state,
            cft,
            model_parameters,
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

function _enzyme_seasonal_joint_loss_block(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
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
    cft = _replace_cft_parameters(base_cft, theta_cft, cft_parameter_names)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            cft,
            model_parameters,
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

function _enzyme_seasonal_soil_loss_block(
    theta_soil::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
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
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            base_cft,
            model_parameters,
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
