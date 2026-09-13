"""
lai_crop!(crop, CFT)

Update leaf-area index from phenology and carbon state.
"""
function lai_crop!(crop,
                   CFT::CFTParameters
)

    launch_1D!(
        lai_crop_kernel!,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_previous_potential,
        crop_prognostic(crop).phenology.senescence,
        crop_prognostic(crop).phenology.senescence_previous,
        crop_prognostic(crop).water.sufficiency,
        crop_prognostic(crop).nitrogen.sufficiency,
        crop_canopy_auxiliary(crop).flaimax,
        crop_prognostic(crop).canopy.laimax_adjusted,
        crop_prognostic(crop).phenology.is_growing,
        CFT,
    )

end

"""
    stress_canopy_loss(stress, rate)

Fraction of standing leaf area a stressed crop sheds today.

Returns exactly zero when `rate` is zero, so the shipped canopy is the inherited
one bit for bit. `stress` is the same `min(wscal, vscal)` that already scales
the expansion increment; this is the other half, the leaf area a stressed crop
loses rather than the leaf area it fails to add.
"""
@inline function stress_canopy_loss(stress::T, rate::T) where {T <: AbstractFloat}
    rate > zero(T) || return zero(T)
    return clamp(rate * (one(T) - clamp(stress, zero(T), one(T))), zero(T), one(T))
end

@kernel inbounds = true function lai_crop_kernel!(
                                  crop_lai::AbstractArray{T},
                                  crop_lai_previous_potential::AbstractArray{T},
                                  crop_senescence::AbstractArray{B},
                                  crop_senescence0::AbstractArray{B},
                                  crop_wscal::AbstractArray{T},
                                  crop_vscal::AbstractArray{T},
                                  crop_flaimax::AbstractArray{T},
                                  crop_laimax_adjusted::AbstractArray{T},
                                  crop_isgrowing::AbstractArray{S},
                                  CFT::CFTParameters
) where {T <: AbstractFloat, S <: Integer, B <: Bool}

    cell = @index(Global)

    @unpack sla, laimax, stress_canopy_loss_rate = CFT

    if crop_isgrowing[cell] == 1
        lai0 = crop_lai[cell]
        stress = min(crop_wscal[cell], crop_vscal[cell])
        retained = one(T) - stress_canopy_loss(stress, T(stress_canopy_loss_rate))
        if !crop_senescence[cell]
            potential_lai = crop_flaimax[cell] * laimax
            # LPJmL's `lai000` stores the previous *potential* LAI, distinct
            # from the actual leaf area retained after water/N limitation.
            # Keeping that state prevents a newly sown winter crop from losing
            # its seed LAI while vernalization keeps potential LAI at zero.
            lai_inc = (potential_lai - crop_lai_previous_potential[cell]) * stress
            crop_lai_previous_potential[cell] = potential_lai
            # The shed fraction leaves the STANDING canopy, not the potential
            # trajectory: `lai_previous_potential` is untouched, so tomorrow's
            # increment is still the phenological one and the crop does not
            # regrow what the drought took.
            crop_lai[cell] = (lai_inc + lai0) * retained
        else
            if !crop_senescence0[cell]
                crop_laimax_adjusted[cell] = crop_lai[cell]
            end
            # Senescence recomputes LAI from `laimax_adjusted` every day, so the
            # loss has to be taken there or it would be undone tomorrow. This is
            # the terminal-drought signature: leaf area stripped faster than the
            # programmed senescence curve.
            crop_laimax_adjusted[cell] *= retained
            crop_lai[cell] = crop_flaimax[cell] * crop_laimax_adjusted[cell]
        end
    else
        crop_lai[cell] = zero(T)
        crop_lai_previous_potential[cell] = zero(T)
        crop_laimax_adjusted[cell] = zero(T)
    end
end
