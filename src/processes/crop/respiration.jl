"""
respiration!(crop, PFT, temp, assim; lpjmlparams=lpjmlparams)

Compute maintenance and growth respiration and update `crop.fluxes.carbon.respiration`.
"""
function respiration_reference!(crop::Crop,
                                PFT::PftParameters,
                                temp::AbstractArray{T},
                                assim::AbstractArray{T};
                                lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    # kernel based
    gtemp_air = crop.workspace.respiration_temperature_response
    launch_1D!(
        temp_response_kernel!,
        temp,
        temp_response,
        e0,
        gtemp_air,
    )

    # unlimited nitrogen
    rosoresp = crop.state.carbon.root * respcoeff * k * nc_ratio.root .* gtemp_air .+ crop.state.carbon.storage * respcoeff * k * nc_ratio.sto .* gtemp_air
    presp = crop.state.carbon.pool * respcoeff * k * nc_ratio.pool .* gtemp_air
    gresp = (assim .- rosoresp .- presp) * r_growth

    # # differentiation based
    # # gate = max.(temp .+ T(40.0), T(0.0)) ./ (max.(temp .+ T(40.0), T(1e-5)))
    # gate = sigmoid.(T(10.0) * (temp .+ T(40.0)))
    # gtemp_air = gate .* exp.(e0 * (one(T) / (temp_response + T(10.0)) .- one(T) ./ (temp .+ temp_response)))
    # rosoresp = crop.state.carbon.root * respcoeff * k .* (crop.state.nitrogen.root ./ (crop.state.carbon.root .+ T(1e-5))) .* gtemp_air .+ crop.state.carbon.storage * respcoeff * k .* (crop.state.nitrogen.storage ./ (crop.state.carbon.storage .+ T(1e-5))) .* gtemp_air
    # presp = crop.state.carbon.pool * respcoeff * k .* (crop.state.nitrogen.pool ./ (crop.state.carbon.pool .+ T(1e-5))) .* gtemp_air
    # gresp = (assim .- rosoresp .- presp) * r_growth
    # gresp = ifelse.(gresp .< zero(T), zero(T), gresp)

    crop.fluxes.carbon.respiration .= (rosoresp .+ presp .+ gresp) .* crop.state.phenology.is_growing

end

# Compatibility/reference entry for callers that already provide net
# assimilation. The daily model uses the allocation-free five-argument path.
function respiration!(crop::Crop,
                      PFT::PftParameters,
                      temp::AbstractArray{T},
                      assim::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    return respiration_reference!(
        crop, PFT, temp, assim; lpjmlparams = lpjmlparams,
    )
end

"""Allocation-free daily respiration using one cell-local CPU/GPU kernel."""
function respiration!(crop::Crop,
                      PFT::PftParameters,
                      temp::AbstractArray{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        respiration_kernel!,
        crop.fluxes.carbon.respiration,
        crop.workspace.respiration_temperature_response,
        crop.state.carbon.root,
        crop.state.carbon.storage,
        crop.state.carbon.pool,
        crop.state.phenology.is_growing,
        temp,
        gross_assimilation,
        leaf_respiration,
        PFT,
        lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function respiration_kernel!(
    respiration::AbstractVector{T},
    temperature_response_buffer::AbstractVector{T},
    root_carbon::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    pool_carbon::AbstractVector{T},
    is_growing::AbstractVector{I},
    temperature::AbstractVector{T},
    gross_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    PFT::PftParameters,
    lpjmlparams::LPJmLParams,
) where {T <: AbstractFloat, I <: Integer}
    cell = @index(Global)
    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    gtemp = temperature[cell] >= T(-40) ? exp(
        T(e0) * (
            one(T) / (T(temp_response) + T(10)) -
            one(T) / (temperature[cell] + T(temp_response))
        ),
    ) : zero(T)
    temperature_response_buffer[cell] = gtemp

    root_storage_respiration =
        root_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.root) * gtemp +
        storage_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.sto) * gtemp
    pool_respiration =
        pool_carbon[cell] * T(respcoeff) * T(k) * T(nc_ratio.pool) * gtemp
    assimilation = gross_assimilation[cell] - leaf_respiration[cell]
    growth_respiration =
        (assimilation - root_storage_respiration - pool_respiration) * T(r_growth)
    respiration[cell] =
        (root_storage_respiration + pool_respiration + growth_respiration) *
        T(is_growing[cell])
end


@kernel inbounds = true function temp_response_kernel!(
                                       temp::AbstractArray{T},
                                       temp_response::T,
                                       e0::T,
                                       gtemp_response::AbstractArray{T}


) where {T <: AbstractFloat}

    cell = @index(Global)

    if temp[cell] >= -40.0
        gtemp_response[cell] = exp(e0 * (one(T) / (temp_response + T(10.0)) - one(T) / (temp[cell] + temp_response)))
    else
        gtemp_response[cell] = zero(T)
    end
end
