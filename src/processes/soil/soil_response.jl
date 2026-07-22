"""
soil_decomp_response!(soil; lpjmlparams=lpjmlparams, soil_decomp_params=soil_decomp_params)

Compute and store shared LPJmL decomposition response terms:
- `soil_decomposition_auxiliary(soil).response` for soil layers;
- row 1 of `litter_response` for surface litter using litter temperature and
  wetness;
- rows 2–3 for incorporated and below-ground litter using the top soil layer.
"""
function soil_decomp_response_reference!(soil;
                                         lpjmlparams::LPJmLParams = lpjmlparams,
                                         soil_decomp_params::SoilDecompParams = soil_decomp_params
)
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params
    moisture = soil_decomposition_workspace(soil).layer_scratch_1
    temperature_response = soil_decomposition_workspace(soil).layer_scratch_2

    # Reuse the two layer scratch buffers for wilting ice, liquid pore
    # capacity, moisture, and temperature response in sequence.
    moisture .= soil_water_auxiliary(soil).wilting_storage .* soil_water_prognostic(soil).wilting_ice_fraction
    temperature_response .= soil_water_auxiliary(soil).saturation_storage .- moisture .-
        soil_water_prognostic(soil).available_ice_storage .- soil_water_prognostic(soil).free_ice_storage
    moisture .= (
        soil_water_auxiliary(soil).relative_content .* soil_water_auxiliary(soil).holding_capacity_storage .+
        soil_water_auxiliary(soil).wilting_storage .- moisture .+ soil_water_auxiliary(soil).free_water
    ) ./ max.(temperature_response, eps)
    moisture .= clamp.(moisture, eps, 1.0f0)

    launch_1D!(
        soil_temperature_response_kernel!,
        soil_thermal_prognostic(soil).temperature,
        lpjmlparams.e0,
        lpjmlparams.temp_response,
        temperature_response,
    )
    soil_decomposition_auxiliary(soil).response .= temperature_response .* (
        intercept .+ moist3 .* moisture.^3 .+
        moist2 .* moisture.^2 .+ moist1 .* moisture
    )
    soil_decomposition_auxiliary(soil).response .= clamp.(soil_decomposition_auxiliary(soil).response, 0.0f0, 1.0f0)

    surface_wetness = soil_decomposition_workspace(soil).surface_scratch_1
    surface_temperature_response = soil_decomposition_workspace(soil).surface_scratch_2
    @views surface_wetness .= ifelse.(
        soil_surface_litter_auxiliary(soil).water_capacity .> eps,
        clamp.(soil_surface_litter_prognostic(soil).water_storage ./
               max.(soil_surface_litter_auxiliary(soil).water_capacity, eps), 0.0f0, 1.0f0),
        moisture[1, :],
    )
    launch_1D!(
        soil_temperature_response_kernel!,
        soil_surface_litter_prognostic(soil).temperature,
        lpjmlparams.e0,
        lpjmlparams.temp_response,
        surface_temperature_response,
    )
    surface_temperature_response .*= (
        intercept .+ moist3 .* surface_wetness.^3 .+
        moist2 .* surface_wetness.^2 .+ moist1 .* surface_wetness
    )
    @views soil_decomposition_auxiliary(soil).litter_response[1, :] .= ifelse.(
        temperature_response[1, :] .> 0.0f0,
        clamp.(surface_temperature_response, 0.0f0, 1.0f0),
        0.0f0,
    )
    @views soil_decomposition_auxiliary(soil).litter_response[2, :] .= soil_decomposition_auxiliary(soil).response[1, :]
    @views soil_decomposition_auxiliary(soil).litter_response[3, :] .= soil_decomposition_auxiliary(soil).response[1, :]
    return nothing
end

function soil_decomp_response!(soil;
                               lpjmlparams::LPJmLParams = lpjmlparams,
                               soil_decomp_params::SoilDecompParams = soil_decomp_params)
    launch_custom!(
        soil_decomp_response_kernel!,
        soil_decomposition_auxiliary(soil).response,
        size(soil_decomposition_auxiliary(soil).response, 2),
        soil_decomposition_auxiliary(soil).litter_response,
        soil_decomposition_workspace(soil).layer_scratch_1,
        soil_decomposition_workspace(soil).layer_scratch_2,
        soil_decomposition_workspace(soil).surface_scratch_1,
        soil_decomposition_workspace(soil).surface_scratch_2,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_prognostic(soil).wilting_ice_fraction,
        soil_water_auxiliary(soil).saturation_storage,
        soil_water_prognostic(soil).available_ice_storage,
        soil_water_prognostic(soil).free_ice_storage,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).holding_capacity_storage,
        soil_water_auxiliary(soil).free_water,
        soil_thermal_prognostic(soil).temperature,
        soil_surface_litter_auxiliary(soil).water_capacity,
        soil_surface_litter_prognostic(soil).water_storage,
        soil_surface_litter_prognostic(soil).temperature,
        eltype(soil_decomposition_auxiliary(soil).response)(lpjmlparams.e0),
        eltype(soil_decomposition_auxiliary(soil).response)(lpjmlparams.temp_response),
        soil_decomp_params,
        size(soil_decomposition_auxiliary(soil).response, 1),
    )
    return nothing
end

@inline function soil_temperature_response_scalar(temperature::T,
                                                  e0::T,
                                                  temp_response::T) where {T <: AbstractFloat}
    if temperature < T(-15)
        return zero(T)
    end
    bounded_temperature = min(temperature, T(40))
    return exp(e0 * (
        one(T) / (temp_response + T(10)) -
        one(T) / (bounded_temperature + temp_response)
    ))
end

@kernel inbounds = true function soil_decomp_response_kernel!(
    response::AbstractMatrix{T},
    litter_response::AbstractMatrix{T},
    moisture_scratch::AbstractMatrix{T},
    temperature_scratch::AbstractMatrix{T},
    surface_wetness_scratch::AbstractVector{T},
    surface_temperature_scratch::AbstractVector{T},
    wilting_storage::AbstractMatrix{T},
    wilting_ice_fraction::AbstractMatrix{T},
    saturation_storage::AbstractMatrix{T},
    available_ice_storage::AbstractMatrix{T},
    free_ice_storage::AbstractMatrix{T},
    relative_content::AbstractMatrix{T},
    holding_capacity_storage::AbstractMatrix{T},
    free_water::AbstractMatrix{T},
    soil_temperature::AbstractMatrix{T},
    surface_water_capacity::AbstractVector{T},
    surface_water_storage::AbstractVector{T},
    surface_temperature::AbstractVector{T},
    e0::T,
    temp_response::T,
    params::SoilDecompParams,
    soil_layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack intercept, moist3, moist2, moist1, eps = params
    top_moisture = zero(T)
    top_temperature_response = zero(T)
    for layer in 1:soil_layers
        wilting_ice = wilting_storage[layer, cell] * wilting_ice_fraction[layer, cell]
        liquid_capacity = saturation_storage[layer, cell] - wilting_ice -
            available_ice_storage[layer, cell] - free_ice_storage[layer, cell]
        moisture = (
            relative_content[layer, cell] * holding_capacity_storage[layer, cell] +
            wilting_storage[layer, cell] - wilting_ice + free_water[layer, cell]
        ) / max(liquid_capacity, T(eps))
        moisture = clamp(moisture, T(eps), one(T))
        temperature_factor = soil_temperature_response_scalar(
            soil_temperature[layer, cell], e0, temp_response,
        )
        moisture_factor = T(intercept) + T(moist3) * moisture^3 +
            T(moist2) * moisture^2 + T(moist1) * moisture
        moisture_scratch[layer, cell] = moisture
        temperature_scratch[layer, cell] = temperature_factor
        response[layer, cell] = clamp(
            temperature_factor * moisture_factor, zero(T), one(T),
        )
        if layer == 1
            top_moisture = moisture
            top_temperature_response = temperature_factor
        end
    end

    surface_wetness = surface_water_capacity[cell] > T(eps) ?
        clamp(surface_water_storage[cell] /
              max(surface_water_capacity[cell], T(eps)), zero(T), one(T)) :
        top_moisture
    surface_temperature_factor = soil_temperature_response_scalar(
        surface_temperature[cell], e0, temp_response,
    )
    surface_moisture_factor = T(intercept) + T(moist3) * surface_wetness^3 +
        T(moist2) * surface_wetness^2 + T(moist1) * surface_wetness
    surface_response = surface_temperature_factor * surface_moisture_factor
    surface_wetness_scratch[cell] = surface_wetness
    surface_temperature_scratch[cell] = surface_response
    litter_response[1, cell] = top_temperature_response > zero(T) ?
        clamp(surface_response, zero(T), one(T)) : zero(T)
    litter_response[2, cell] = response[1, cell]
    litter_response[3, cell] = response[1, cell]
end

@kernel inbounds = true function soil_temperature_response_kernel!(
    temperature::AbstractArray{T},
    e0::T,
    temp_response::T,
    response::AbstractArray{T},
) where {T <: AbstractFloat}
    index = @index(Global)
    if temperature[index] >= T(-15)
        bounded_temperature = min(temperature[index], T(40))
        response[index] = exp(e0 * (
            one(T) / (temp_response + T(10)) -
            one(T) / (bounded_temperature + temp_response)
        ))
    else
        response[index] = zero(T)
    end
end


"""
temp_response(temp; lpjmlparams=lpjmlparams)

Calculate the temperature response function for soil decomposition based on
LPJmL's formulation. The response is zero below `-15 °C`, uses the
Lloyd–Taylor expression from `-15 °C` through `40 °C`, and is capped at
its `40 °C` value above that threshold.
"""
function temp_response(temp::AbstractArray{T};
                       lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack e0, temp_response = lpjmlparams

    bounded_temp = clamp.(temp, T(-15.0), T(40.0))
    response = exp.(e0 .* (
        one(T) / (temp_response + T(10.0)) .-
        one(T) ./ (bounded_temp .+ temp_response)
    ))
    return ifelse.(temp .>= T(-15.0), response, zero(T))
end
