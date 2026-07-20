"""
crop_carbon!(photos, PFT, crop, pet, soil, temp, co2)

Run the daily crop carbon process chain: respiration, allocation, and phenology coupling.
"""
function crop_carbon!(photos::CropPhotosynthesis,
                      crop::Crop,
                      output::Output,
                      PFT::PftParameters,
                      temp::AbstractArray{T};
                      output_row::Union{Nothing, Integer} = nothing
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(
        crop, PFT, temp,
        photos.gross_assimilation,
        photos.leaf_respiration,
    )

    # compute crop carbon allocation
    carbon_allocation!(PFT, crop, photos)
    # crop.carbon.organs = vcat(reshape(crop.carbon.root, (1, :)), reshape(crop.carbon.leaf, (1, :)), reshape(crop.carbon.storage, (1, :)), reshape(crop.carbon.pool, (1, :)))

    sources = (
        gpp = photos.gross_assimilation,
        npp = crop.carbon.npp,
        lambda = photos.lambda,
        potential_vmax = photos.potential_vmax,
        vmax = photos.vmax,
        nitrogen_limitation = photos.nitrogen_limitation,
        respiration = crop.carbon.respiration,
        lai = crop.canopy.lai,
        fphu = crop.phenology.fphu,
        biomass = crop.carbon.biomass,
    )
    for (field, source) in pairs(sources)
        if output_row === nothing
            setproperty!(
                output.crop,
                field,
                _append_output_row(getproperty(output.crop, field), source),
            )
        else
            _write_output_row!(getproperty(output.crop, field), output_row, source)
        end
    end

end
