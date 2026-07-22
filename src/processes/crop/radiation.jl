"""
petpar!(pet, day, lat, temp, lwnet, swdown; dayseconds=86400)

Compute daylength, PAR, and equilibrium evapotranspiration diagnostics.
"""
function petpar_reference!(pet::PetPar,
                           day::Int64,
                           lat::AbstractArray{T},
                           temp::AbstractArray{T},
                           lwnet::AbstractArray{T},
                           swdown::AbstractArray{T};
                           dayseconds = 86400
) where {T <: AbstractFloat}


    delta = T(deg2rad(-23.4 * cos(2 * π * (day + 10) / 365)))
    u = T.(sin.(deg2rad.(lat)) * sin(delta))
    v = T.(cos.(deg2rad.(lat)) * cos(delta))

    launch_1D!(
        daylength_kernel!,
        pet.daylength,
        u,
        v,
    )

    swnet = (1 .- pet.albedo) .* swdown

    pet.par .= T(dayseconds) .* swdown ./ T(2)

    s = T(2.503e6) * exp.(T(17.269) * temp ./ (T(237.3) .+ temp)) ./
        ((T(237.3) .+ temp).^2)

    gamma_t = T(65.05) .+ T(0.064) * temp
    lambda = T(2.495e6) .- T(2380) * temp

    pet.eeq .= T(dayseconds) * (s ./ (s .+ gamma_t) ./ lambda) .*
        (swnet .+ lwnet .* (pet.daylength / T(24)))

    # idx = pet.eeq .< 0
    # pet.eeq[idx] .= zero(T)
    # pet.eeq .= ifelse.(pet.eeq .< 0, zero(T), pet.eeq)
    pet.eeq .= max.(pet.eeq, zero(T))

    ## check equilibrium evapotranspiration
    pet.eeq .= min.(pet.eeq, T(15)) ##  set an upper bound for pet.eeq to avoid extreme values to stop GPU computing

end

"""Allocation-free radiation and equilibrium-evaporation preprocessing."""
function petpar!(pet::PetPar,
                 day::Int64,
                 lat::AbstractArray{T},
                 temp::AbstractArray{T},
                 lwnet::AbstractArray{T},
                 swdown::AbstractArray{T};
                 dayseconds = 86400
) where {T <: AbstractFloat}
    delta = T(deg2rad(-23.4 * cos(2 * π * (day + 10) / 365)))
    launch_1D!(
        petpar_kernel!,
        pet.daylength,
        pet.par,
        pet.eeq,
        pet.albedo,
        lat,
        temp,
        lwnet,
        swdown,
        delta,
        T(dayseconds),
    )
    return nothing
end

@kernel inbounds = true function petpar_kernel!(
    daylength::AbstractVector{T},
    par::AbstractVector{T},
    eeq::AbstractVector{T},
    albedo::AbstractVector{T},
    latitude::AbstractVector{T},
    temperature::AbstractVector{T},
    lwnet::AbstractVector{T},
    swdown::AbstractVector{T},
    delta::T,
    dayseconds::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    # Retain the original vector path's Float64 angular conversion before
    # narrowing u/v to the model precision. This keeps daylength numerically
    # identical while all large state arrays remain Float32.
    latitude_radians = Float64(latitude[cell]) * π / 180.0
    u = T(sin(latitude_radians) * sin(delta))
    v = T(cos(latitude_radians) * cos(delta))
    daylight = if u >= v
        T(24)
    elseif u <= -v
        zero(T)
    else
        T(T(24) * acos(-u / v) * (1 / π))
    end
    daylength[cell] = daylight
    par[cell] = dayseconds * swdown[cell] / T(2)

    temperature_offset = T(237.3) + temperature[cell]
    slope = T(2.503e6) * exp(
        T(17.269) * temperature[cell] / temperature_offset,
    ) / (temperature_offset^2)
    psychrometric = T(65.05) + T(0.064) * temperature[cell]
    latent_heat = T(2.495e6) - T(2380) * temperature[cell]
    net_shortwave = (one(T) - albedo[cell]) * swdown[cell]
    equilibrium = dayseconds * slope / (slope + psychrometric) / latent_heat *
        (net_shortwave + lwnet[cell] * daylight / T(24))
    eeq[cell] = clamp(equilibrium, zero(T), T(15))
end

@kernel inbounds = true function daylength_kernel!(
                                   pet_daylength::AbstractArray{T},
                                   u::AbstractArray{T},
                                   v::AbstractArray{T}
) where {T <: AbstractFloat}

    cell = @index(Global)

    if u[cell] >= v[cell]
        pet_daylength[cell] = 24
    elseif u[cell] <= -v[cell]
        pet_daylength[cell] = 0
    else
        hh = acos(-u[cell] / v[cell])
        pet_daylength[cell] = 24 * hh * (1 / π)
    end
end

# for one cft
"""
apar_crop!(PFT, crop, pet)

Compute absorbed PAR and fPAR for non-maize crops.
"""
function apar_crop_reference!(PFT::PftParameters,
                              crop,
                              pet::PetPar,
                              snow_height = nothing,
)

    @unpack name, lightextcoeff, albedo_leaf, alphaa  = PFT

    actual_lai = max.(
        zero(eltype(crop_prognostic(crop).canopy.lai)),
        crop_prognostic(crop).canopy.lai .- crop_prognostic(crop).canopy.lai_npp_deficit,
    )
    crop_canopy_auxiliary(crop).fpar .= 1 .- exp.(-lightextcoeff * actual_lai)
    if snow_height !== nothing
        crop_canopy_auxiliary(crop).fpar .*= snow_height .<= zero(eltype(snow_height))
    end

    crop_canopy_auxiliary(crop).apar .= pet.par * (1 - albedo_leaf) * alphaa .* crop_canopy_auxiliary(crop).fpar

end

function apar_crop!(PFT::PftParameters, crop, pet::PetPar, snow_height = nothing)
    T = eltype(crop_canopy_auxiliary(crop).apar)
    launch_1D!(
        apar_crop_kernel!,
        crop_canopy_auxiliary(crop).apar,
        crop_canopy_auxiliary(crop).fpar,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        snow_height === nothing ? pet.eeq : snow_height,
        pet.par,
        T(PFT.lightextcoeff),
        T(PFT.albedo_leaf),
        T(PFT.alphaa),
        false,
        snow_height !== nothing,
    )
    return nothing
end
# Radiation and daylength preprocessing for canopy photosynthesis.

"""
apar_crop_maize!(PFT, crop, pet)

Compute absorbed PAR and maize-specific fPAR parameterization.
"""
function apar_crop_maize_reference!(PFT::PftParameters,
                                    crop,
                                    pet::PetPar,
                                    snow_height = nothing,
)

    @unpack name, lightextcoeff, albedo_leaf, alphaa  = PFT

    actual_lai = max.(
        zero(eltype(crop_prognostic(crop).canopy.lai)),
        crop_prognostic(crop).canopy.lai .- crop_prognostic(crop).canopy.lai_npp_deficit,
    )
    crop_canopy_auxiliary(crop).fpar .= min.(1.0f0, max.(0.0f0, 0.2558f0 * max.(0.01f0, actual_lai) .- 0.0024f0))
    if snow_height !== nothing
        crop_canopy_auxiliary(crop).fpar .*= snow_height .<= zero(eltype(snow_height))
    end

    crop_canopy_auxiliary(crop).apar .= pet.par * (1 - albedo_leaf) * alphaa .* crop_canopy_auxiliary(crop).fpar

end


function apar_crop_maize!(PFT::PftParameters, crop, pet::PetPar, snow_height = nothing)
    T = eltype(crop_canopy_auxiliary(crop).apar)
    launch_1D!(
        apar_crop_kernel!,
        crop_canopy_auxiliary(crop).apar,
        crop_canopy_auxiliary(crop).fpar,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        snow_height === nothing ? pet.eeq : snow_height,
        pet.par,
        T(PFT.lightextcoeff),
        T(PFT.albedo_leaf),
        T(PFT.alphaa),
        true,
        snow_height !== nothing,
    )
    return nothing
end

@kernel inbounds = true function apar_crop_kernel!(
    apar::AbstractVector{T},
    fpar::AbstractVector{T},
    lai::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    snow_height::AbstractVector{T},
    par::AbstractVector{T},
    light_extinction::T,
    leaf_albedo::T,
    alpha_a::T,
    maize::Bool,
    apply_snow_cover::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    actual_lai = max(zero(T), lai[cell] - lai_npp_deficit[cell])
    absorbed_fraction = if maize
        min(one(T), max(zero(T), T(0.2558) * max(T(0.01), actual_lai) - T(0.0024)))
    else
        one(T) - exp(-light_extinction * actual_lai)
    end
    snow_free = !apply_snow_cover || snow_height[cell] <= zero(T)
    absorbed_fraction *= T(snow_free)
    fpar[cell] = absorbed_fraction
    apar[cell] = par[cell] * (one(T) - leaf_albedo) * alpha_a * absorbed_fraction
end
