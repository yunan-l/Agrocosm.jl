"""
crop_carbon!(photos, PFT, crop, pet, soil, temp, co2)

Run the daily crop carbon process chain: respiration, allocation, and phenology coupling.
"""
function crop_carbon!(photos::CropPhotosynthesis,
                      crop::Crop,
                      output::Output,
                      PFT::PftParameters,
                      temp::AbstractArray{T}
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(crop, PFT, temp, photos.gross_assimilation - photos.leaf_respiration)

    # compute crop carbon allocation
    carbon_allocation!(PFT, crop, photos)
    # crop.carbon.organs = vcat(reshape(crop.carbon.root, (1, :)), reshape(crop.carbon.leaf, (1, :)), reshape(crop.carbon.storage, (1, :)), reshape(crop.carbon.pool, (1, :)))

    output.crop.gpp = vcat(output.crop.gpp, reshape(photos.gross_assimilation, (1, :)))
    output.crop.npp = vcat(output.crop.npp, reshape(crop.carbon.npp, (1, :)))
    output.crop.lai = vcat(output.crop.lai, reshape(crop.canopy.lai, (1, :)))
    output.crop.fphu = vcat(output.crop.fphu, reshape(crop.phenology.fphu, (1, :)))
    output.crop.biomass = vcat(output.crop.biomass, reshape(crop.carbon.biomass, (1, :)))

end