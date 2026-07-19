
"""
albedo!(PFT, crop, pet)

Update canopy albedo terms used by PET and radiation partitioning.
"""
function albedo!(PFT::PftParameters,
                 crop::Crop,
                 pet::PetPar;
                 soil_albedo = 0.3f0  # Albedo of bare soil (0-1). Should be soil and soil moisture dependent */
)

    @unpack fpc = PFT

    crop_albedo!(PFT, crop)

    pet.albedo .= crop.canopy.albedo .+ (1 .- fpc * crop.phenology.is_growing) * soil_albedo

end


function crop_albedo!(PFT::PftParameters,
                      crop::Crop
)
    @unpack albedo_leaf, albedo_litter, fpc = PFT

    albedo_green_leaves = fpc * crop.canopy.phenology_fraction * albedo_leaf

    # albedo of PFT without green foliage (litter background albedo)

    albedo_brown_litter = fpc * (1 .- crop.canopy.phenology_fraction) * albedo_litter

    crop.canopy.albedo .= albedo_green_leaves .+ albedo_brown_litter

end
