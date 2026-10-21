@generated function _replace_cft_parameters(
    base_cft::Agrocosm.CFTParameters{T, S},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, S, parameter_names}
    fields = fieldnames(Agrocosm.CFTParameters)
    values = map(fields) do field
        if field === :temp_co2
            index = findfirst(==(:temp_co2_high), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :temp_co2)) :
                :(Agrocosm.TempCO2{$T}(
                    getfield(getfield(base_cft, :temp_co2), :low),
                    theta[$index],
                ))
        elseif field === :temp_photos
            index = findfirst(==(:temp_photos_high), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :temp_photos)) :
                :(Agrocosm.TempPhotos{$T}(
                    getfield(getfield(base_cft, :temp_photos), :low),
                    theta[$index],
                ))
        elseif field === :basetemp
            index = findfirst(==(:basetemp_low), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :basetemp)) :
                :(Agrocosm.BaseTemp{$T}(
                    theta[$index],
                    getfield(getfield(base_cft, :basetemp), :high),
                ))
        end
        index = findfirst(==(field), parameter_names)
        index === nothing ? :(getfield(base_cft, $(QuoteNode(field)))) : :(theta[$index])
    end
    return :(Agrocosm.CFTParameters{$T, $S}($(values...)))
end

@generated function _replace_cft_shadow(
    zero_cft::Agrocosm.CFTParameters{T, S},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, S, parameter_names}
    fields = fieldnames(Agrocosm.CFTParameters)
    values = map(fields) do field
        if field === :temp_co2
            index = findfirst(==(:temp_co2_high), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :temp_co2)) :
                :(Agrocosm.TempCO2{$T}(
                    getfield(getfield(zero_cft, :temp_co2), :low),
                    dtheta[$index],
                ))
        elseif field === :temp_photos
            index = findfirst(==(:temp_photos_high), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :temp_photos)) :
                :(Agrocosm.TempPhotos{$T}(
                    getfield(getfield(zero_cft, :temp_photos), :low),
                    dtheta[$index],
                ))
        elseif field === :basetemp
            index = findfirst(==(:basetemp_low), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :basetemp)) :
                :(Agrocosm.BaseTemp{$T}(
                    dtheta[$index],
                    getfield(getfield(zero_cft, :basetemp), :high),
                ))
        end
        index = findfirst(==(field), parameter_names)
        index === nothing ? :(getfield(zero_cft, $(QuoteNode(field)))) : :(dtheta[$index])
    end
    return :(Agrocosm.CFTParameters{$T, $S}($(values...)))
end

_replace_cft_shadow(zero_cft, dtheta, parameter_names::Tuple) =
    _replace_cft_shadow(zero_cft, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_cft_parameters)},
    ::Type,
    base_cft::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_cft_parameters(base_cft.val, theta.val, parameter_names.val)
    shadow = _replace_cft_shadow(
        Enzyme.make_zero(base_cft.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_cft_parameters(
    base_cft::Agrocosm.CFTParameters{T, S},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::Agrocosm.CFTParameters{T, S} where {T, S}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_cft_parameters(base_cft, theta, Val(parameter_names))
end

"""Replace only the selected soil-related fields in `LPJmLParams`."""
@generated function _replace_lpjml_parameters(
    base_lpjml::Agrocosm.LPJmLParams{T},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    fields = fieldnames(Agrocosm.LPJmLParams)
    values = map(fields) do field
        if field === :k_soil10
            index = findfirst(==(:k_soil10_fast), parameter_names)
            fast = index === nothing ?
                :(getfield(getfield(base_lpjml, :k_soil10), :fast)) :
                :(theta[$index])
            return :(Agrocosm.K_Soil10{$T}(
                $fast,
                getfield(getfield(base_lpjml, :k_soil10), :slow),
            ))
        end
        index = findfirst(==(field), parameter_names)
        return index === nothing ?
            :(getfield(base_lpjml, $(QuoteNode(field)))) :
            :(theta[$index])
    end
    return :(Agrocosm.LPJmLParams{$T}($(values...)))
end

@generated function _replace_lpjml_shadow(
    zero_lpjml::Agrocosm.LPJmLParams{T},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    fields = fieldnames(Agrocosm.LPJmLParams)
    values = map(fields) do field
        if field === :k_soil10
            index = findfirst(==(:k_soil10_fast), parameter_names)
            fast = index === nothing ?
                :(getfield(getfield(zero_lpjml, :k_soil10), :fast)) :
                :(dtheta[$index])
            return :(Agrocosm.K_Soil10{$T}(
                $fast,
                getfield(getfield(zero_lpjml, :k_soil10), :slow),
            ))
        end
        index = findfirst(==(field), parameter_names)
        return index === nothing ?
            :(getfield(zero_lpjml, $(QuoteNode(field)))) :
            :(dtheta[$index])
    end
    return :(Agrocosm.LPJmLParams{$T}($(values...)))
end

_replace_lpjml_shadow(zero_lpjml, dtheta, parameter_names::Tuple) =
    _replace_lpjml_shadow(zero_lpjml, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_lpjml_parameters)},
    ::Type,
    base_lpjml::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_lpjml_parameters(
        base_lpjml.val, theta.val, parameter_names.val,
    )
    shadow = _replace_lpjml_shadow(
        Enzyme.make_zero(base_lpjml.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_lpjml_parameters(
    base_lpjml::Agrocosm.LPJmLParams{T},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::Agrocosm.LPJmLParams{T} where {T}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_lpjml_parameters(base_lpjml, theta, Val(parameter_names))
end

"""Rebuild ModelParameters while keeping all non-soil bundles fixed."""
@generated function _replace_model_parameters(
    base_parameters::Agrocosm.ModelParameters{T},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    return :(Agrocosm.ModelParameters{$T}(
        _replace_lpjml_parameters(
            getfield(base_parameters, :lpjml),
            theta,
            Val($(QuoteNode(parameter_names))),
        ),
        getfield(base_parameters, :photosynthesis),
        getfield(base_parameters, :snow),
        getfield(base_parameters, :soil_thermal),
        getfield(base_parameters, :soil_decomposition),
    ))
end

@generated function _replace_model_shadow(
    zero_parameters::Agrocosm.ModelParameters{T},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    return :(Agrocosm.ModelParameters{$T}(
        _replace_lpjml_shadow(
            getfield(zero_parameters, :lpjml),
            dtheta,
            Val($(QuoteNode(parameter_names))),
        ),
        getfield(zero_parameters, :photosynthesis),
        getfield(zero_parameters, :snow),
        getfield(zero_parameters, :soil_thermal),
        getfield(zero_parameters, :soil_decomposition),
    ))
end

_replace_model_shadow(zero_parameters, dtheta, parameter_names::Tuple) =
    _replace_model_shadow(zero_parameters, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_model_parameters)},
    ::Type,
    base_parameters::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_model_parameters(
        base_parameters.val, theta.val, parameter_names.val,
    )
    shadow = _replace_model_shadow(
        Enzyme.make_zero(base_parameters.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_model_parameters(
    base_parameters::Agrocosm.ModelParameters{T},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::Agrocosm.ModelParameters{T} where {T}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_model_parameters(base_parameters, theta, Val(parameter_names))
end
