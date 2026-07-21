"""
crop_carbon!(photos, PFT, crop, pet, soil, temp, co2)

Run the daily crop carbon process chain: respiration, allocation, and phenology coupling.
"""
function crop_carbon!(crop::Crop,
                      output::Output,
                      PFT::PftParameters,
                      air_temperature::AbstractVector{T},
                      soil_temperature::AbstractMatrix{T};
                      output_row::Union{Nothing, Integer} = nothing,
                      lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(
        crop, PFT, air_temperature, soil_temperature,
        crop.fluxes.carbon.gross_assimilation,
        crop.fluxes.carbon.leaf_respiration;
        lpjmlparams = lpjmlparams,
    )

    # compute crop carbon allocation
    carbon_allocation!(PFT, crop)

    sources = (
        gpp = crop.fluxes.carbon.gross_assimilation,
        npp = crop.fluxes.carbon.npp,
        lambda = crop.auxiliary.photosynthesis.lambda,
        potential_vcmax = crop.auxiliary.photosynthesis.potential_vcmax,
        vcmax = crop.auxiliary.photosynthesis.vcmax,
        nitrogen_limitation = crop.auxiliary.photosynthesis.nitrogen_limitation,
        respiration = crop.fluxes.carbon.respiration,
        lai = crop.auxiliary.canopy.actual_lai,
        fphu = crop.auxiliary.phenology.fphu,
        biomass = crop.state.carbon.biomass,
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
