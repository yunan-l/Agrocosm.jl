"""
respiration!(crop, PFT, temp, assim; lpjmlparams=lpjmlparams)

Compute maintenance and growth respiration and update `crop.fluxes.carbon.respiration`.
"""
function respiration_reference!(crop::Crop,
                                PFT::PftParameters,
                                air_temperature::AbstractVector{T},
                                soil_temperature::AbstractMatrix{T},
                                assim::AbstractArray{T};
                                lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    bounded_air = min.(air_temperature, T(40))
    bounded_soil = min.(vec(@view(soil_temperature[1, :])), T(40))
    gtemp_air = ifelse.(
        air_temperature .>= T(-15),
        exp.(T(e0) .* (one(T) / (T(temp_response) + T(10)) .-
                       one(T) ./ (bounded_air .+ T(temp_response)))),
        zero(T),
    )
    gtemp_soil = ifelse.(
        vec(@view(soil_temperature[1, :])) .>= T(-15),
        exp.(T(e0) .* (one(T) / (T(temp_response) + T(10)) .-
                       one(T) ./ (bounded_soil .+ T(temp_response)))),
        zero(T),
    )
    # unlimited nitrogen
    roresp = crop.state.carbon.root * respcoeff * k * nc_ratio.root .* gtemp_soil
    soresp = crop.state.carbon.storage * respcoeff * k * nc_ratio.sto .* gtemp_air
    presp = crop.state.carbon.pool * respcoeff * k * nc_ratio.pool .* gtemp_air
    gresp = max.(zero(T), (assim .- roresp .- soresp .- presp) * r_growth)

    # # differentiation based
    # # gate = max.(temp .+ T(40.0), T(0.0)) ./ (max.(temp .+ T(40.0), T(1e-5)))
    # gate = sigmoid.(T(10.0) * (temp .+ T(40.0)))
    # gtemp_air = gate .* exp.(e0 * (one(T) / (temp_response + T(10.0)) .- one(T) ./ (temp .+ temp_response)))
    # rosoresp = crop.state.carbon.root * respcoeff * k .* (crop.state.nitrogen.root ./ (crop.state.carbon.root .+ T(1e-5))) .* gtemp_air .+ crop.state.carbon.storage * respcoeff * k .* (crop.state.nitrogen.storage ./ (crop.state.carbon.storage .+ T(1e-5))) .* gtemp_air
    # presp = crop.state.carbon.pool * respcoeff * k .* (crop.state.nitrogen.pool ./ (crop.state.carbon.pool .+ T(1e-5))) .* gtemp_air
    # gresp = (assim .- rosoresp .- presp) * r_growth
    # gresp = ifelse.(gresp .< zero(T), zero(T), gresp)

    crop.fluxes.carbon.respiration .=
        (roresp .+ soresp .+ presp .+ gresp) .* crop.state.phenology.is_growing
end

# Compatibility/reference entry for callers that already provide net
# assimilation. The daily model uses the allocation-free five-argument path.
function respiration!(crop::Crop,
                      PFT::PftParameters,
                      temp::AbstractArray{T},
                      assim::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    return respiration_reference!(crop, PFT, temp, reshape(temp, 1, :), assim;
                                  lpjmlparams = lpjmlparams)
end

"""Allocation-free daily respiration using one cell-local CPU/GPU kernel."""
function respiration!(crop::Crop,
                      PFT::PftParameters,
                      air_temperature::AbstractVector{T},
                      soil_temperature::AbstractMatrix{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        respiration_kernel!,
        crop.fluxes.carbon.respiration,
        crop.state.carbon.root,
        crop.state.carbon.storage,
        crop.state.carbon.pool,
        crop.state.phenology.is_growing,
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
function respiration!(crop::Crop,
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
