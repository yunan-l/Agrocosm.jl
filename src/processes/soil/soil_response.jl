"""
soil_decomp_response!(soil; lpjmlparams=lpjmlparams, soil_decomp_params=soil_decomp_params)

Compute and store shared soil decomposition response terms:
- `soil.decomposition.response` for soil layers
- `soil.decomposition.litter_response` for litter decomposition (using top layer response)
"""
function soil_decomp_response!(soil::Soil;
                               lpjmlparams::LPJmLParams = lpjmlparams,
                               soil_decomp_params::SoilDecompParams = soil_decomp_params
)
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    moist = (soil.water.relative_content .* soil.water.holding_capacity_storage .+ soil.water.wilting_storage .+ soil.water.free_water) ./ max.(soil.water.saturation_storage, eps) # soil.water.saturation_storage .- soil.water.wilting_storage
    moist = clamp.(moist, eps, 1.0f0)
    gtemp_soil = temp_response(soil.thermal.temperature; lpjmlparams = lpjmlparams)

    soil.decomposition.response .= gtemp_soil .* (intercept .+ moist3 .* moist.^3 .+ moist2 .* moist.^2 .+ moist1 .* moist)
    soil.decomposition.response .= clamp.(soil.decomposition.response, 0.0f0, 1.0f0)

    soil.decomposition.litter_response .= reshape(soil.decomposition.response[1, :], (1, :))
end


"""
temp_response(temp; lpjmlparams=lpjmlparams)

Calculate the temperature response function for soil decomposition based on LPJmL's formulation. The function is defined as
`g(T) = exp(e0 * (1/(temp_response+10) - 1/(T+temp_response)))` for `T >= -40`,
otherwise `0`.
"""
function temp_response(temp::AbstractArray{T};
                       lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack e0, temp_response = lpjmlparams

    return ifelse.(temp .>= T(-40.0), exp.(e0 .* (one(T) / (temp_response + T(10.0)) .- one(T) ./ (temp .+ temp_response))), zero(T))
end


