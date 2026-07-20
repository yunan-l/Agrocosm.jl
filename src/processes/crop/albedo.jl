
"""
albedo!(PFT, crop, pet)

Update canopy albedo terms used by PET and radiation partitioning.
"""
function albedo_reference!(PFT::PftParameters,
                           crop::Crop,
                           pet::PetPar;
                           soil_albedo = 0.3f0
)

    @unpack fpc = PFT

    crop_albedo!(PFT, crop)

    pet.albedo .= crop.auxiliary.canopy.albedo .+ (1 .- fpc .* crop.state.phenology.is_growing) .* soil_albedo

end

function albedo!(PFT::PftParameters,
                 crop::Crop,
                 pet::PetPar;
                 soil_albedo = 0.3f0)
    T = eltype(pet.albedo)
    launch_1D!(
        albedo_kernel!,
        pet.albedo,
        crop.auxiliary.canopy.albedo,
        crop.auxiliary.canopy.phenology_fraction,
        crop.state.phenology.is_growing,
        T(PFT.albedo_leaf),
        T(PFT.albedo_litter),
        T(PFT.fpc),
        T(soil_albedo),
    )
    return nothing
end

@kernel inbounds = true function albedo_kernel!(
    pet_albedo::AbstractVector{T},
    canopy_albedo::AbstractVector{T},
    phenology_fraction::AbstractVector{T},
    is_growing::AbstractVector{S},
    leaf_albedo::T,
    litter_albedo::T,
    fpc::T,
    soil_albedo::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    canopy = fpc * (
        phenology_fraction[cell] * leaf_albedo +
        (one(T) - phenology_fraction[cell]) * litter_albedo
    )
    canopy_albedo[cell] = canopy
    pet_albedo[cell] = canopy +
        (one(T) - fpc * is_growing[cell]) * soil_albedo
end


function crop_albedo!(PFT::PftParameters,
                      crop::Crop
)
    @unpack albedo_leaf, albedo_litter, fpc = PFT

    # Fuse green-leaf and brown-litter terms into the preallocated canopy buffer.
    crop.auxiliary.canopy.albedo .= fpc .* (
        crop.auxiliary.canopy.phenology_fraction .* albedo_leaf .+
        (1 .- crop.auxiliary.canopy.phenology_fraction) .* albedo_litter
    )

end
