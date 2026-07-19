"""
respiration!(crop, PFT, temp, assim; lpjmlparams=lpjmlparams)

Compute maintenance and growth respiration and update `crop.carbon.respiration`.
"""
function respiration!(crop::Crop,
                      PFT::PftParameters,
                      temp::AbstractArray{T},
                      assim::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    # kernel based
    gtemp_air = similar(temp)
    launch_1D!(
        temp_response_kernel!,
        temp,
        temp_response,
        e0,
        gtemp_air,
    )

    # unlimited nitrogen
    rosoresp = crop.carbon.root * respcoeff * k * nc_ratio.root .* gtemp_air .+ crop.carbon.storage * respcoeff * k * nc_ratio.sto .* gtemp_air
    presp = crop.carbon.pool * respcoeff * k * nc_ratio.pool .* gtemp_air
    gresp = (assim .- rosoresp .- presp) * r_growth

    # # differentiation based
    # # gate = max.(temp .+ T(40.0), T(0.0)) ./ (max.(temp .+ T(40.0), T(1e-5)))
    # gate = sigmoid.(T(10.0) * (temp .+ T(40.0)))
    # gtemp_air = gate .* exp.(e0 * (one(T) / (temp_response + T(10.0)) .- one(T) ./ (temp .+ temp_response)))
    # rosoresp = crop.carbon.root * respcoeff * k .* (crop.nitrogen.root ./ (crop.carbon.root .+ T(1e-5))) .* gtemp_air .+ crop.carbon.storage * respcoeff * k .* (crop.nitrogen.storage ./ (crop.carbon.storage .+ T(1e-5))) .* gtemp_air
    # presp = crop.carbon.pool * respcoeff * k .* (crop.nitrogen.pool ./ (crop.carbon.pool .+ T(1e-5))) .* gtemp_air
    # gresp = (assim .- rosoresp .- presp) * r_growth
    # gresp = ifelse.(gresp .< zero(T), zero(T), gresp)

    crop.carbon.respiration = (rosoresp .+ presp .+ gresp) .* crop.phenology.is_growing

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
