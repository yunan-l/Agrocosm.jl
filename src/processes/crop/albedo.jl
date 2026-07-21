
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
        crop.state.canopy.lai,
        crop.state.phenology.is_growing,
        T(PFT.laimax),
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
    lai::AbstractVector{T},
    is_growing::AbstractVector{S},
    laimax::T,
    leaf_albedo::T,
    litter_albedo::T,
    fpc::T,
    soil_albedo::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    phenology_fraction = lai[cell] / laimax
    canopy = fpc * (
        phenology_fraction * leaf_albedo +
        (one(T) - phenology_fraction) * litter_albedo
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
    phenology_fraction = crop.state.canopy.lai ./ PFT.laimax
    crop.auxiliary.canopy.albedo .= fpc .* (
        phenology_fraction .* albedo_leaf .+
        (1 .- phenology_fraction) .* albedo_litter
    )

end
