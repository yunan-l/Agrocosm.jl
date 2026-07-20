"""
soil_decomp_response!(soil; lpjmlparams=lpjmlparams, soil_decomp_params=soil_decomp_params)

Compute and store shared LPJmL decomposition response terms:
- `soil.decomposition.response` for soil layers;
- row 1 of `litter_response` for surface litter using litter temperature and
  wetness;
- rows 2–3 for incorporated and below-ground litter using the top soil layer.
"""
function soil_decomp_response!(soil::Soil;
                               lpjmlparams::LPJmLParams = lpjmlparams,
                               soil_decomp_params::SoilDecompParams = soil_decomp_params
)
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    wilting_ice = soil.water.wilting_storage .* soil.water.wilting_ice_fraction
    liquid_pwp = soil.water.wilting_storage .- wilting_ice
    liquid_pore_capacity = soil.water.saturation_storage .-
        wilting_ice .-
        soil.water.available_ice_storage .-
        soil.water.free_ice_storage
    moist = (soil.water.relative_content .* soil.water.holding_capacity_storage .+
             liquid_pwp .+ soil.water.free_water) ./
            max.(liquid_pore_capacity, eps)
    moist = clamp.(moist, eps, 1.0f0)
    gtemp_soil = temp_response(soil.thermal.temperature; lpjmlparams = lpjmlparams)

    soil.decomposition.response .= gtemp_soil .* (intercept .+ moist3 .* moist.^3 .+ moist2 .* moist.^2 .+ moist1 .* moist)
    soil.decomposition.response .= clamp.(soil.decomposition.response, 0.0f0, 1.0f0)

    surface_wetness = ifelse.(
        soil.surface_litter.water_capacity .> eps,
        clamp.(soil.surface_litter.water_storage ./
               max.(soil.surface_litter.water_capacity, eps), 0.0f0, 1.0f0),
        vec(@view moist[1, :]),
    )
    surface_temperature_response = temp_response(
        soil.surface_litter.temperature;
        lpjmlparams = lpjmlparams,
    )
    surface_response = surface_temperature_response .* (
        intercept .+ moist3 .* surface_wetness.^3 .+
        moist2 .* surface_wetness.^2 .+ moist1 .* surface_wetness
    )
    soil.decomposition.litter_response[1, :] .=
        clamp.(surface_response, 0.0f0, 1.0f0)
    soil.decomposition.litter_response[2, :] .=
        @view soil.decomposition.response[1, :]
    soil.decomposition.litter_response[3, :] .=
        @view soil.decomposition.response[1, :]
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
