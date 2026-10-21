"""Compute a checkpointed reverse gradient for the independent soil block."""
function Agrocosm.enzyme_seasonal_soil_gradient_blockwise(
    theta_soil::AbstractVector{T},
    state_factory::F,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat, F}
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))

    state, _ = _state_and_shadow(state_factory)
    Agrocosm.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_soil_loss_block(
            theta_soil,
            state,
            base_cft,
            base_parameters,
            soil_parameter_names,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            nitrogen_limit_vcmax,
            crop_resp_fix,
        )
    end

    gradient = zeros(T, length(theta_soil))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            Agrocosm.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = zeros(T, length(theta_soil))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_soil_loss_block,
            Enzyme.Duplicated(theta_soil, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(base_parameters),
            Enzyme.Const(soil_parameter_names),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax),
            Enzyme.Const(crop_resp_fix),
        )
        reverse_primal += result[2]
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end

    return (
        primal = reverse_primal,
        forward_primal = forward_primal,
        gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

"""Compute one checkpointed reverse gradient for CFT and soil parameters."""
function Agrocosm.enzyme_seasonal_joint_gradient_blockwise(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state_factory::F,
    base_cft::Agrocosm.CFTParameters,
    base_parameters::Agrocosm.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat, F}
    length(theta_cft) == length(cft_parameter_names) || throw(DimensionMismatch(
        "CFT parameter values and names must have identical lengths",
    ))
    length(theta_soil) == length(soil_parameter_names) || throw(DimensionMismatch(
        "soil parameter values and names must have identical lengths",
    ))
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))

    state, _ = _state_and_shadow(state_factory)
    Agrocosm.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_joint_loss_block(
            theta_cft,
            theta_soil,
            state,
            base_cft,
            base_parameters,
            cft_parameter_names,
            soil_parameter_names,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            nitrogen_limit_vcmax,
            crop_resp_fix,
        )
    end

    cft_gradient = zeros(T, length(theta_cft))
    soil_gradient = zeros(T, length(theta_soil))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            Agrocosm.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_cft_gradient = zeros(T, length(theta_cft))
        block_soil_gradient = zeros(T, length(theta_soil))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_joint_loss_block,
            Enzyme.Duplicated(theta_cft, block_cft_gradient),
            Enzyme.Duplicated(theta_soil, block_soil_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(base_parameters),
            Enzyme.Const(cft_parameter_names),
            Enzyme.Const(soil_parameter_names),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax),
            Enzyme.Const(crop_resp_fix),
        )
        reverse_primal += result[2]
        cft_gradient .+= block_cft_gradient
        soil_gradient .+= block_soil_gradient
        state_cotangent = block_shadow
    end

    return (;
        primal = reverse_primal,
        forward_primal,
        cft_gradient,
        soil_gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

function _seasonal_block_ranges(days::Tuple, block_days::Int)
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    ranges = UnitRange{Int}[]
    first_day = 1
    while first_day <= length(days)
        last_day = min(first_day + block_days - 1, length(days))
        push!(ranges, first_day:last_day)
        first_day = last_day + 1
    end
    return ranges
end

"""
    Agrocosm.enzyme_seasonal_gradient_blockwise(
        theta, state_factory, base_cft, global_parameters, climate,
        parameter_names, days, layer_depth, context; block_days=30,
    )

Compute the exact reverse gradient of a fixed seasonal loss with explicit
periodic state checkpoints. A forward pass stores one state snapshot at each
block boundary. Blocks are then differentiated from the final block backwards;
the state shadow from a later block seeds the output-state cotangent of the
preceding block. This is intentionally separate from the production runtime.
"""
function Agrocosm.enzyme_seasonal_gradient_blockwise(
    theta::AbstractVector{T},
    state_factory::F,
    base_cft::Agrocosm.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::Agrocosm.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat, F}
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))

    state, _ = _state_and_shadow(state_factory)
    Agrocosm.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_loss_block(
            theta,
            state,
            base_cft,
            global_parameters,
            climate,
            parameter_names,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            nitrogen_limit_vcmax,
            crop_resp_fix,
        )
    end

    gradient = zeros(T, length(theta))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            Agrocosm.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = zeros(T, length(theta))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_loss_block,
            Enzyme.Duplicated(theta, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(global_parameters),
            Enzyme.Const(climate),
            Enzyme.Const(parameter_names),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax),
            Enzyme.Const(crop_resp_fix),
        )
        reverse_primal += result[2]
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end

    return (
        primal = reverse_primal,
        forward_primal = forward_primal,
        gradient = gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

"""
    Agrocosm.enzyme_daily_transition_objective(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, day, observable, layer_depth,
    )

Run the opt-in continuous portion of one existing C3 daily transition and
return one scalar daily observable. The caller should prepare the lifecycle
state with discrete events already fixed; output diagnostics are intentionally
excluded from this AD path.
"""
function Agrocosm.enzyme_daily_transition_objective(
    theta::AbstractVector{T},
    state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    day::Integer,
    observable::Symbol,
    layer_depth;
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    crop_resp_fix::Bool = true,
) where {T <: AbstractFloat}
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    _enzyme_continuous_transition!(
        state, cft, global_parameters, climate, day, observable, layer_depth,
        irrigation, nitrogen_limit_vcmax, crop_resp_fix, nitrogen_limit_vcmax,
    )
    return _daily_transition_observable(state, observable)
end
