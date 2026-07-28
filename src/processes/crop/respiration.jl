"""
respiration!(crop, PFT, temp, assim; lpjmlparams=lpjmlparams)

Compute maintenance and growth respiration and update `crop_fluxes(crop).carbon.respiration`.
"""

"""Allocation-free daily respiration using one cell-local CPU/GPU kernel."""
function respiration!(crop,
                      PFT::PftParameters,
                      air_temperature::AbstractVector{T},
                      soil_temperature::AbstractMatrix{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        respiration_kernel!,
        crop_fluxes(crop).carbon.respiration,
        crop_prognostic(crop).carbon.root,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.pool,
        crop_prognostic(crop).phenology.is_growing,
        air_temperature,
        soil_temperature,
        gross_assimilation,
        leaf_respiration,
        PFT,
        lpjmlparams,
    )
    return nothing
end

# Backward-compatible entry for callers without an explicit soil-temperature
# profile. Daily simulations use the method above.
function respiration!(crop,
                      PFT::PftParameters,
                      air_temperature::AbstractVector{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams) where {T <: AbstractFloat}
    return respiration!(
        crop, PFT, air_temperature, reshape(air_temperature, 1, :),
        gross_assimilation, leaf_respiration; lpjmlparams = lpjmlparams,
    )
end

@kernel inbounds = true function respiration_kernel!(
    respiration::AbstractVector{T},
    root_carbon::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    pool_carbon::AbstractVector{T},
    is_growing::AbstractVector{I},
    air_temperature::AbstractVector{T},
    soil_temperature::AbstractMatrix{T},
    gross_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    PFT::PftParameters,
    lpjmlparams::LPJmLParams,
) where {T <: AbstractFloat, I <: Integer}
    cell = @index(Global)
    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    air_temp = min(air_temperature[cell], T(40))
    soil_temp = min(soil_temperature[1, cell], T(40))
    gtemp_air = air_temperature[cell] >= T(-15) ? exp(
        T(e0) * (one(T) / (T(temp_response) + T(10)) -
                 one(T) / (air_temp + T(temp_response))),
    ) : zero(T)
    gtemp_soil = soil_temperature[1, cell] >= T(-15) ? exp(
        T(e0) * (one(T) / (T(temp_response) + T(10)) -
                 one(T) / (soil_temp + T(temp_response))),
    ) : zero(T)
    root_respiration =
        root_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.root) * gtemp_soil
    storage_respiration =
        storage_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.sto) * gtemp_air
    pool_respiration =
        pool_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.pool) * gtemp_air
    assimilation = gross_assimilation[cell] - leaf_respiration[cell]
    growth_respiration = max(
        zero(T),
        (assimilation - root_respiration - storage_respiration - pool_respiration) *
        T(r_growth),
    )
    active = T(is_growing[cell])
    respiration[cell] =
        (root_respiration + storage_respiration + pool_respiration + growth_respiration) * active
end
