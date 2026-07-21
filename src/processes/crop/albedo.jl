"""
    albedo!(PFT, crop, soil, pet; maize=false)

Compute LPJmL-style effective surface albedo from green canopy absorption,
surface-litter cover, bare soil, and the snow state present at the start of the
day. Surface-litter cover is rebuilt directly from its carbon stock so event
routing and restart initialization cannot leave radiation on a stale cache.
"""
function albedo_reference!(PFT::PftParameters,
                           crop::Crop,
                           soil::Soil,
                           pet::PetPar;
                           maize::Bool = false,
                           soil_albedo = 0.3f0,
                           snow_albedo = 0.65f0,
                           litter_carbon_fraction = 0.42f0)
    T = eltype(pet.albedo)
    actual_lai = max.(
        zero(T),
        crop.state.canopy.lai .- crop.state.canopy.lai_npp_deficit,
    )
    green_fraction = if maize
        clamp.(T(0.2558) .* max.(T(0.01), actual_lai) .- T(0.0024), zero(T), one(T))
    else
        one(T) .- exp.(-T(PFT.lightextcoeff) .* actual_lai)
    end
    litter_dry_matter = max.(@view(soil.carbon.litter[1, :]), zero(T)) ./
        T(litter_carbon_fraction)
    litter_cover = one(T) .- exp.(-T(6e-3) .* litter_dry_matter)
    snow_cover = soil.snow.height .> zero(T)

    green_albedo = green_fraction .* ifelse.(
        snow_cover, T(snow_albedo), T(PFT.albedo_leaf),
    )
    background_fraction = one(T) .- green_fraction
    litter_albedo = litter_cover .* background_fraction .* ifelse.(
        snow_cover, T(snow_albedo), T(PFT.albedo_litter),
    )
    soil_background = (one(T) .- litter_cover) .* background_fraction
    soil_component = ifelse.(
        snow_cover,
        soil_background .* soil.snow.fraction .* T(snow_albedo),
        soil_background .* T(soil_albedo),
    )
    crop_surface_albedo = green_albedo .+ litter_albedo .+ soil_component
    bare_surface_albedo = soil.snow.fraction .* T(snow_albedo) .+
        (one(T) .- soil.snow.fraction) .* T(soil_albedo)
    growing = crop.state.phenology.is_growing .!= 0

    crop.auxiliary.canopy.albedo .= ifelse.(
        growing, crop_surface_albedo, zero(T),
    )
    pet.albedo .= ifelse.(
        growing,
        crop_surface_albedo .+ max(one(T) - T(PFT.fpc), zero(T)) .* bare_surface_albedo,
        bare_surface_albedo,
    )
    return nothing
end

function albedo!(PFT::PftParameters,
                 crop::Crop,
                 soil::Soil,
                 pet::PetPar;
                 maize::Bool = false,
                 soil_albedo = 0.3f0,
                 snow_albedo = 0.65f0,
                 litter_carbon_fraction = 0.42f0)
    T = eltype(pet.albedo)
    launch_1D!(
        albedo_kernel!,
        pet.albedo,
        crop.auxiliary.canopy.albedo,
        crop.state.canopy.lai,
        crop.state.canopy.lai_npp_deficit,
        crop.state.phenology.is_growing,
        soil.carbon.litter,
        soil.snow.height,
        soil.snow.fraction,
        T(PFT.lightextcoeff),
        T(PFT.albedo_leaf),
        T(PFT.albedo_litter),
        T(PFT.fpc),
        T(soil_albedo),
        T(snow_albedo),
        T(litter_carbon_fraction),
        maize,
    )
    return nothing
end

@kernel inbounds = true function albedo_kernel!(
    pet_albedo::AbstractVector{T},
    canopy_albedo::AbstractVector{T},
    lai::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    is_growing::AbstractVector{S},
    carbon_litter::AbstractMatrix{T},
    snow_height::AbstractVector{T},
    snow_fraction::AbstractVector{T},
    light_extinction::T,
    leaf_albedo::T,
    litter_albedo::T,
    fpc::T,
    soil_albedo::T,
    snow_albedo::T,
    litter_carbon_fraction::T,
    maize::Bool,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    actual_lai = max(zero(T), lai[cell] - lai_npp_deficit[cell])
    green_fraction = maize ?
        clamp(T(0.2558) * max(T(0.01), actual_lai) - T(0.0024), zero(T), one(T)) :
        one(T) - exp(-light_extinction * actual_lai)
    litter_dry_matter = max(carbon_litter[1, cell], zero(T)) /
        litter_carbon_fraction
    litter_cover = one(T) - exp(-T(6e-3) * litter_dry_matter)
    background_fraction = one(T) - green_fraction
    snow_present = snow_height[cell] > zero(T)

    green_component = green_fraction *
        (snow_present ? snow_albedo : leaf_albedo)
    litter_component = litter_cover * background_fraction *
        (snow_present ? snow_albedo : litter_albedo)
    soil_background = (one(T) - litter_cover) * background_fraction
    soil_component = snow_present ?
        soil_background * snow_fraction[cell] * snow_albedo :
        soil_background * soil_albedo
    crop_surface = green_component + litter_component + soil_component
    bare_surface = snow_fraction[cell] * snow_albedo +
        (one(T) - snow_fraction[cell]) * soil_albedo

    if is_growing[cell] != zero(S)
        canopy_albedo[cell] = crop_surface
        pet_albedo[cell] = crop_surface + max(one(T) - fpc, zero(T)) * bare_surface
    else
        canopy_albedo[cell] = zero(T)
        pet_albedo[cell] = bare_surface
    end
end
