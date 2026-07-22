"""
lai_crop!(crop, PFT)

Update leaf-area index from phenology and carbon state.
"""
function lai_crop!(crop,
                   PFT::PftParameters
)

    launch_1D!(
        lai_crop_kernel!,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).phenology.senescence,
        crop_prognostic(crop).phenology.senescence_previous,
        crop_prognostic(crop).water.sufficiency,
        crop_prognostic(crop).nitrogen.sufficiency,
        crop_canopy_auxiliary(crop).flaimax,
        crop_prognostic(crop).canopy.laimax_adjusted,
        crop_prognostic(crop).phenology.is_growing,
        PFT,
    )

end

@kernel inbounds = true function lai_crop_kernel!(
                                  crop_lai::AbstractArray{T},
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
    else
        crop_lai[cell] = zero(T)
        crop_laimax_adjusted[cell] = zero(T)
    end
end
