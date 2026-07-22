"""
transpiration!(photos_adtmm, PFT, crop, pet, soil, co2; lpjmlparams=lpjmlparams)

Compute water demand/supply balance and layer-resolved transpiration uptake.
"""
function transpiration!(photos_adtmm::AbstractArray{T},
                        PFT::PftParameters,
                        crop,
                        pet::PetPar,
                        soil,
                        co2::AbstractArray{T};
                        lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    # Root-zone weighted water availability is accumulated inside the cell
    # kernel, avoiding a separate broadcast and reduction array every day.
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop_prognostic(crop).carbon.root))
    # demand = ifelse.(crop_canopy_auxiliary(crop).canopy_conductance .> 0, (1 .- crop_canopy_auxiliary(crop).canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop_canopy_auxiliary(crop).canopy_conductance), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(water_demand_supply_kernel!,
               crop_canopy_auxiliary(crop).canopy_conductance,
               photos_adtmm,
               co2,
               pet.daylength,
               crop_canopy_auxiliary(crop).fpar,
               crop_fluxes(crop).water.transpiration_layer,
               crop_prognostic(crop).water.demand_sum,
               crop_prognostic(crop).water.supply_sum,
               crop_stress_auxiliary(crop).water_deficit,
               crop_prognostic(crop).water.sufficiency,
               crop_prognostic(crop).carbon.root,
               crop_canopy_auxiliary(crop).canopy_wet,
               crop_prognostic(crop).phenology.is_growing,
               pet.eeq,
               crop_root_input(crop).distribution,
               crop_root_auxiliary(crop).zone_available_water,
               soil_water_auxiliary(soil).relative_content,
               soil_water_auxiliary(soil).holding_capacity_storage,
               PFT,
               kernel_params)

end

@kernel inbounds = true function water_demand_supply_kernel!(
                                             crop_gp::AbstractArray{T},
                                             photos_adtmm::AbstractArray{T},
                                             co2::AbstractArray{T},
                                             daylength::AbstractArray{T},
                                             crop_fpar::AbstractArray{T},
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
                                             crop_rootzone_available_water::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             PFT::PftParameters,
                                             kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers = kernel_params

    @unpack ALPHAM, GM, LAMBDA_OPT = lpjmlparams
    @unpack fpc, emax, gmin = PFT

    co2_index = length(co2) == 1 ? 1 : cell
    co2_bar = co2[co2_index] * T(1e-5)
    if co2_bar > zero(T) && daylength[cell] > zero(T)
        conductance_denominator = co2_bar * (one(T) - T(LAMBDA_OPT)) *
            hour2sec(daylength[cell])
        crop_gp[cell] = T(1.6) * photos_adtmm[cell] / conductance_denominator +
            T(gmin) * crop_fpar[cell]
    else
        crop_gp[cell] = zero(T)
    end

    wr = zero(T)
    rootzone_water = zero(T)
    for l in 1:soil_layers
        wr += soil_w[l, cell] * crop_rootdist[l]
        if l <= 3
            rootzone_water += soil_w[l, cell] * soil_whcs[l, cell] * crop_rootdist[l]
        end
    end
    crop_rootzone_available_water[cell] = rootzone_water

    if crop_isgrowing[cell] == 1
        supply = emax * wr * (1 - exp(-0.04f0 * crop_rootc[cell]))
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
            crop_wdf[cell] = clamp(
                T(100.0) * crop_w_supplysum[cell] / crop_w_demandsum[cell],
                zero(T), T(100.0),
            )
        else
            crop_wdf[cell] = T(100.0)
        end

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            crop_wscal[cell] = (emax * wr) / (pet_eeq[cell] * ALPHAM / (one(T) + (GM * ALPHAM) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        # Potential transpiration constrained by demand/supply and canopy fraction.
        if wr > 0
            transp = min(supply, demand) / wr * fpc
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

        if wr > 0
            transp = transp_cor / wr
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
        # Neutral stress for an absent stand; is_growing still gates all fluxes.
        crop_wscal[cell] = one(T)
    end
end
