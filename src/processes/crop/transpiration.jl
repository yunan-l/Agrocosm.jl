"""
transpiration!(photos_adtmm, PFT, crop, pet, soil, co2; lpjmlparams=lpjmlparams)

Compute water demand/supply balance and layer-resolved transpiration uptake.
"""
function transpiration!(photos_adtmm::AbstractArray{T},
                        PFT::PftParameters,
                        crop::Crop,
                        pet::PetPar,
                        soil::Soil,
                        co2::AbstractArray{T};
                        lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack LAMBDA_OPT = lpjmlparams
    @unpack gmin = PFT

    # `co2` is already a partial pressure in Pa. Since 1 Pa = 1e-5 bar,
    # this is LPJmL's ppm2bar(original_co2) without applying ppm conversion twice.
    co2_bar = co2 .* T(1e-5)
    conductance_denominator = co2_bar .* (one(T) - T(LAMBDA_OPT)) .* hour2sec(pet.daylength)
    crop.gp .= ifelse.(
        (co2_bar .> zero(T)) .& (pet.daylength .> zero(T)),
        T(1.6) .* photos_adtmm ./ conductance_denominator .+ T(gmin) .* crop.fpar,
        zero(T),
    )

    # Root-zone weighted soil water availability per cell.
    wr = sum(soil.w .* crop.rootdist, dims = 1)
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop.rootc))
    # demand = ifelse.(crop.gp .> 0, (1 .- crop.canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop.gp), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(water_demand_supply_kernel!,
               crop.gp,
               crop.trans_layer,
               crop.w_demandsum,
               crop.w_supplysum,
               crop.wdf,
               crop.wscal,
               crop.rootc,
               crop.canopy_wet,
               crop.isgrowing,
               pet.eeq,
               crop.rootdist,
               soil.w,
               soil.whcs,
               wr,
               PFT,
               kernel_params)

end

@kernel inbounds = true function water_demand_supply_kernel!(
                                             crop_gp::AbstractArray{T},
                                             crop_trans_layer::AbstractArray{T},
                                             crop_w_demandsum::AbstractArray{T},
                                             crop_w_supplysum::AbstractArray{T},
                                             crop_wdf::AbstractArray{T},
                                             crop_wscal::AbstractArray{T},
                                             crop_rootc::AbstractArray{T},
                                             crop_canopy_wet::AbstractArray{T},
                                             crop_isgrowing::AbstractArray{S},
                                             pet_eeq::AbstractArray{T},
                                             crop_rootdist::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             wr::AbstractArray{T},
                                             PFT::PftParameters,
                                             kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack lpjmlparams, soil_layers = kernel_params

    @unpack ALPHAM, GM = lpjmlparams
    @unpack fpc, emax = PFT

    if crop_isgrowing[cell] == 1
        supply = emax * wr[cell] * (1 - exp(-0.04f0 * crop_rootc[cell]))
        if crop_gp[cell] > 0
            demand = (1 - crop_canopy_wet[cell]) * pet_eeq[cell] * ALPHAM / (1 + (GM * ALPHAM) / crop_gp[cell])
        else
            demand = zero(T)
        end
        
        crop_w_demandsum[cell] += demand
        if supply > demand
            crop_w_supplysum[cell] += demand
        else
            crop_w_supplysum[cell] += supply
        end

        if crop_w_demandsum[cell] > 0.0
            crop_wdf[cell] = T(100.0) * crop_w_supplysum[cell] / crop_w_demandsum[cell]
        else
            crop_wdf[cell] = T(100.0)
        end

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            crop_wscal[cell] = (emax * wr[cell]) / (pet_eeq[cell] * ALPHAM / (one(T) + (GM * ALPHAM) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        # Potential transpiration constrained by demand/supply and canopy fraction.
        if wr[cell] > 0
            transp = min(supply, demand) / wr[cell] * fpc
        else
            transp = zero(T)
        end

        transp_cor = zero(T)

        # Apply layer-wise extraction cap so uptake does not exceed layer storage.
        if transp > 0
            for l in 1:soil_layers
                transp_frac = 1
                if transp * crop_rootdist[l] * soil_w[l, cell] > soil_w[l, cell] * soil_whcs[l, cell]
                    transp_frac = soil_whcs[l, cell] / (transp * crop_rootdist[l])
                end
                transp_tmp = transp * crop_rootdist[l] * soil_w[l, cell] * transp_frac
                if transp_tmp > soil_w[l, cell] * soil_whcs[l, cell]
                    transp_cor += soil_w[l, cell] * soil_whcs[l, cell]
                    if transp_cor < 1.0f-5
                        transp_cor = zero(T)
                    end
                else
                    transp_cor += transp_tmp
                end
            end
        else
            transp_cor = zero(T)
        end

        if wr[cell] > 0
            transp = transp_cor / wr[cell]
        else
            transp = zero(T)
        end

        # LPJmL recomputes actual canopy conductance after layer extraction.
        # Store it in `gp`; downstream lambda solving consumes this actual value.
        actual_supply = fpc > zero(T) ? transp_cor / fpc : zero(T)
        if actual_supply < demand && pet_eeq[cell] > zero(T)
            denominator = (one(T) - crop_canopy_wet[cell]) * pet_eeq[cell] * ALPHAM - actual_supply
            crop_gp[cell] = denominator > zero(T) ?
                            (GM * ALPHAM) * actual_supply / denominator : zero(T)
        end

        # Distribute corrected transpiration back to layers by root distribution.
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = transp * crop_rootdist[l] * soil_w[l, cell]
            if crop_trans_layer[l, cell] > soil_w[l, cell] * soil_whcs[l, cell]
                crop_trans_layer[l, cell] = soil_w[l, cell] * soil_whcs[l, cell]
            end
        end
    else
        crop_gp[cell] = zero(T)
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = zero(T)
        end
        crop_w_demandsum[cell] = zero(T)
        crop_w_supplysum[cell] = zero(T)
        crop_wdf[cell] = zero(T)
        crop_wscal[cell] = zero(T)
    end
end
