# Output buffers are diagnostics for this scalar transition objective. They
# are still executed by the primal function, but do not participate in its
# derivative. Keeping these rules in the optional extension avoids changing
# the production output path or its ordinary runtime behavior.
function Enzyme.EnzymeRules.inactive(
    ::typeof(Agrocosm.prepare_output_block!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(Agrocosm._write_output_row!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(Agrocosm.record_ecosystem_flux_outputs!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(Agrocosm.convert_precision),
    ::Type{<:AbstractFloat},
    ::Agrocosm.ModelParameters,
)
    return nothing
end

# Climate forcing is a fixed input to every current Enzyme objective. Execute
# the production reader for its primal weather-buffer updates, while keeping
# its backend launch outside Enzyme's differentiated graph.
function Enzyme.EnzymeRules.inactive(::typeof(Agrocosm.readclimate!), args...)
    return nothing
end

# Deposition is a fixed additive input to the mineral-N state for the current
# objectives. The primal kernel runs normally; its parameter tangent is the
# unchanged incoming mineral-N tangent.
function Enzyme.EnzymeRules.inactive(::typeof(Agrocosm.nitrogen_deposition!), args...)
    return nothing
end

function _zero_tangent_arrays!(value)
    value isa AbstractArray && return fill!(value, zero(eltype(value)))
    value === nothing && return nothing
    value isa Union{Symbol, String, Function, Type} && return value
    isstructtype(typeof(value)) || return value
    for field in fieldnames(typeof(value))
        _zero_tangent_arrays!(getfield(value, field))
    end
    return value
end

"""Create an Enzyme zero tangent for a mutable Agrocosm state container."""
function Agrocosm.enzyme_zero_tangent(value)
    # Some nested immutable containers receive shallow array fields from
    # `make_zero`; isolate them before clearing the tangent arrays.
    tangent = deepcopy(Enzyme.make_zero(value))
    _zero_tangent_arrays!(tangent)
    return tangent
end

function _state_and_shadow(state_factory)
    pair = state_factory()
    pair isa Tuple && length(pair) == 2 || throw(ArgumentError(
        "state_factory must return (state, state_shadow)",
    ))
    state, state_shadow = pair
    typeof(state) === typeof(state_shadow) || throw(ArgumentError(
        "state and state_shadow must have identical types",
    ))
    return state, state_shadow
end

"""
    Agrocosm.enzyme_forward_directional(
        objective, theta, direction, state_factory, args...;
        runtime_activity=true,
    )

Differentiate a scalar, in-place Agrocosm objective in one parameter direction.
`objective(theta, state, args...)` must return a scalar and may mutate `state`.
`state_factory()` is called for every probe and must return a fresh
`(state, state_shadow)` pair. This keeps the AD call from reusing a mutated
state between parameter directions.
"""
function Agrocosm.enzyme_forward_directional(
    objective::F,
    theta::AbstractVector{T},
    direction::AbstractVector{T},
    state_factory::SF,
    args...;
    runtime_activity::Bool = true,
    return_primal::Bool = false,
) where {F, T <: AbstractFloat, SF}
    axes(theta) == axes(direction) || throw(DimensionMismatch(
        "theta and direction must have identical axes",
    ))
    state, state_shadow = _state_and_shadow(state_factory)
    mode = if return_primal
        Enzyme.set_runtime_activity(Enzyme.ForwardWithPrimal)
    elseif runtime_activity
        Enzyme.set_runtime_activity(Enzyme.Forward)
    else
        Enzyme.Forward
    end
    annotations = map(Enzyme.Const, args)
    result = Enzyme.autodiff(
        mode,
        objective,
        Enzyme.Duplicated,
        Enzyme.Duplicated(theta, direction),
        Enzyme.Duplicated(state, state_shadow),
        annotations...,
    )
    return return_primal ? (primal = result[2], directional = result[1]) : result[1]
end

"""
    Agrocosm.enzyme_forward_gradient(objective, theta, state_factory, args...)

Compute one forward-mode Enzyme derivative for each element of `theta`. The
state factory is called once per parameter, so each derivative starts from the
same state. This is intentionally a small validation API; reverse-mode
optimization will be added only after the in-place state contract is verified.
"""
function Agrocosm.enzyme_forward_gradient(
    objective::F,
    theta::AbstractVector{T},
    state_factory::SF,
    args...;
    runtime_activity::Bool = true,
) where {F, T <: AbstractFloat, SF}
    gradient = similar(theta)
    direction = similar(theta)
    for index in eachindex(theta)
        fill!(direction, zero(T))
        direction[index] = one(T)
        gradient[index] = Agrocosm.enzyme_forward_directional(
            objective,
            theta,
            direction,
            state_factory,
            args...;
            runtime_activity,
        )
    end
    return gradient
end

"""Preallocate the one-row output buffer required by the AD daily wrapper."""
function Agrocosm.enzyme_prepare_daily_state!(state::Agrocosm.ModelState)
    Agrocosm.prepare_output_block!(state.output, 1, 1; reuse = true)
    return state
end

# `crop_carbon!` combines continuous carbon processes with dynamic diagnostic
# field writes. Reuse its two process calls while leaving the diagnostic path
# outside the Enzyme graph.
function _enzyme_crop_carbon!(
    state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters,
    air_temperature,
    soil_temperature,
    lpjmlparams,
    include_biological_fixation_cost::Bool,
    crop_resp_fix::Bool,
)
    fluxes = Agrocosm.crop_fluxes(state)
    Agrocosm.respiration!(
        state,
        cft,
        air_temperature,
        soil_temperature,
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration;
        crop_resp_fix,
        include_biological_fixation_cost,
        lpjmlparams,
    )
    Agrocosm.carbon_allocation!(
        cft, state; include_biological_fixation_cost,
    )
    return nothing
end
