"""
lai_crop!(crop, PFT)

Update leaf-area index from phenology and carbon state.
"""
function lai_crop!(crop::Crop,
                   PFT::PftParameters
)

    launch_1D!(
        lai_crop_kernel!,
        crop.canopy.lai,
        crop.canopy.phenology_fraction,
        crop.phenology.senescence,
        crop.phenology.senescence_previous,
        crop.water.stress,
        crop.nitrogen.stress,
        crop.canopy.flaimax,
        crop.canopy.laimax_adjusted,
        crop.phenology.is_growing,
        PFT,
    )

end

@kernel inbounds = true function lai_crop_kernel!(
                                  crop_lai::AbstractArray{T},
                                  phenology_fraction::AbstractArray{T},
                                  crop_senescence::AbstractArray{B},
                                  crop_senescence0::AbstractArray{B},
                                  crop_wscal::AbstractArray{T},
                                  crop_vscal::AbstractArray{T},
                                  crop_flaimax::AbstractArray{T},
                                  crop_laimax_adjusted::AbstractArray{T},
                                  crop_isgrowing::AbstractArray{S},
                                  PFT::PftParameters
) where {T <: AbstractFloat, S <: Integer, B <: Bool}

    cell = @index(Global)

    @unpack sla, laimax = PFT

    if crop_isgrowing[cell] == 1
        lai0 = crop_lai[cell]
        if !crop_senescence[cell]
            crop_lai[cell] = crop_flaimax[cell] * laimax
            # scale daily LAI increment with minimum of wscal and vscal as simplest approach
            lai_inc = (crop_lai[cell] - lai0) * min(crop_wscal[cell]/T(1.5), crop_vscal[cell])
            crop_lai[cell] = lai_inc + lai0
        else
            if !crop_senescence0[cell]
                crop_laimax_adjusted[cell] = crop_lai[cell]
            end
            crop_lai[cell] = crop_flaimax[cell] * crop_laimax_adjusted[cell]
        end
        # if !crop_senescence[cell]
        #     # scale daily LAI increment with minimum of wscal and vscal as simplest approach
        #     lai_inc = (crop_lai[cell] - lai0) * min(crop_wscal[cell]/T(1.5), crop_vscal[cell])
        #     crop_lai[cell] = lai_inc + lai0
        # end
    else
        crop_lai[cell] = zero(T)
        crop_laimax_adjusted[cell] = zero(T)
    end
    phenology_fraction[cell] = crop_lai[cell] / laimax
end


"""
lai_deficit!(crop, PFT)

Apply LAI deficit correction under senescence or carbon-limited states.
"""
function lai_deficit!(crop::Crop,
                      PFT::PftParameters
)

    launch_1D!(
        lai_deficit_kernel!,
        crop.canopy.lai,
        crop.phenology.senescence,
        crop.carbon.biomass,
        crop.carbon.root,
        crop.carbon.leaf,
        crop.canopy.lai_npp_deficit,
        crop.phenology.is_growing,
        PFT,
    )

end

@kernel inbounds = true function lai_deficit_kernel!(
                                     crop_lai::AbstractArray{T},
                                     crop_senescence::AbstractArray{B},
                                     crop_biomass::AbstractArray{T},
                                     crop_rootc::AbstractArray{T},
                                     crop_leafc::AbstractArray{T},
                                     crop_lai_nppdeficit::AbstractArray{T},
                                     crop_isgrowing::AbstractArray{S},
                                     PFT::PftParameters
) where {T <: AbstractFloat, S <: Integer, B <: Bool}

    cell = @index(Global)

    @unpack sla = PFT

    if crop_isgrowing[cell] == 1
        if !crop_senescence[cell]
            if (crop_biomass[cell] - crop_rootc[cell]) >= crop_lai[cell] / sla
                crop_lai_nppdeficit[cell] = zero(T)
            else
                crop_lai_nppdeficit[cell] = crop_lai[cell] - crop_leafc[cell] * sla
                # today's lai_deficit is subtracted from tomorrow's LAI in lai_crop(), fpar_crop(), and actual_lai_crop().
                # These routines account for LAI effects on the simulation.
            end
        end
    else
        crop_lai_nppdeficit[cell] = zero(T)
    end
end
