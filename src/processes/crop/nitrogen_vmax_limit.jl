"""
    limit_vmax_by_nitrogen!(crop, PFT, temperature)

Apply LPJmL's crop leaf-nitrogen constraint to the potential Rubisco capacity.
`crop.nitrogen.demand_leaf` is the leaf-N stock remaining after uptake logic;
structural leaf N is protected and only excess N supports Rubisco activity.
The result is bounded to `[eps(T), potential_vmax]` for an active crop and the
dimensionless retained fraction is stored in `photos.nitrogen_limitation`.
"""
function limit_vmax_by_nitrogen!(crop::Crop,
                                 PFT::PftParameters,
                                 temperature::AbstractArray{T};
                                 lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        nitrogen_vmax_limit_kernel!,
        crop.photosynthesis.vmax,
        crop.photosynthesis.potential_vmax,
        crop.photosynthesis.nitrogen_limitation,
        crop.nitrogen.demand_leaf,
        crop.carbon.leaf,
        crop.phenology.is_growing,
        temperature,
        PFT,
        lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function nitrogen_vmax_limit_kernel!(
    vmax::AbstractArray{T},
    potential_vmax::AbstractArray{T},
    nitrogen_limitation::AbstractArray{T},
    available_leaf_nitrogen::AbstractArray{T},
    leaf_carbon::AbstractArray{T},
    is_growing::AbstractArray{S},
    temperature::AbstractArray{T},
    PFT::PftParameters,
    lpjmlparams::LPJmLParams,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    potential = max(zero(T), potential_vmax[cell])

    if is_growing[cell] == one(S) && potential > zero(T)
        @unpack ncleaf = PFT
        @unpack p, k_temp = lpjmlparams
        rubisco_nitrogen = max(
            zero(T),
            available_leaf_nitrogen[cell] - T(ncleaf.low) * leaf_carbon[cell],
        )
        nitrogen_capacity = rubisco_nitrogen /
            exp(-T(k_temp) * (temperature[cell] - T(25))) /
            (T(p) * T(1e-3)) * (T(86400) * T(12) * T(1e-6))
        limited = min(potential, max(eps(T), nitrogen_capacity))
        vmax[cell] = limited
        nitrogen_limitation[cell] = clamp(limited / potential, zero(T), one(T))
    else
        vmax[cell] = potential
        nitrogen_limitation[cell] = potential > zero(T) ? one(T) : zero(T)
    end
end
