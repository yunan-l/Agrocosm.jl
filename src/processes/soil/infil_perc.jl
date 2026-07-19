"""
infil_perc!(soil; lpjmlparams=lpjmlparams)

Update infiltration, percolation, runoff, and nitrate transport through soil layers.
"""
function infil_perc!(soil::Soil;
                    lpjmlparams::LPJmLParams = lpjmlparams
)
    # One-cell kernel launch; each thread updates the full vertical soil column for that cell.
    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5, anion_excl = 0.3f0, NPERCO = 0.4f0)

    launch_1D!(infil_perc_kernel!,
               soil.water.infiltration,
               soil.water.relative_content,
               soil.water.holding_capacity_storage,
               soil.water.free_water,
               soil.water.ice_storage,
               soil.water.available_ice_storage,
               soil.water.free_ice_storage,
               soil.water.saturation_fraction,
               soil.water.saturation_storage,
               soil.water.wilting_fraction,
               soil.water.wilting_storage,
               soil.water.influx,
               soil.water.outflux,
               soil.water.saturated_conductivity,
               soil.water.surface_runoff,
               soil.water.lateral_runoff,
               soil.water.bottom_drainage,
               soil.water.percolation,
               soil.thermal.freeze_depth,
               soil.surface_litter.cover,
               soil.water.beta,
               soil.nitrogen.nitrate,
               soil.nitrogen.leaching,
               soil.properties.layer_depth,
               kernel_params)

end

@kernel inbounds = true function infil_perc_kernel!(
                                    infil::AbstractArray{T},
                                    soil_w::AbstractArray{M},
                                    soil_whcs::AbstractArray{M},
                                    soil_w_fw::AbstractArray{M},
                                    soil_ice::AbstractArray{M},
                                    soil_ice_available::AbstractArray{M},
                                    soil_ice_free::AbstractArray{M},
                                    soil_wsat::AbstractArray{M},
                                    soil_wsats::AbstractArray{M},
                                    soil_wpwp::AbstractArray{M},
                                    soil_wpwps::AbstractArray{M},
                                    soil_w_influx::AbstractArray{M},
                                    soil_w_outflux::AbstractArray{M},
                                    soil_Ks::AbstractArray{M},
                                    soil_srunoff::AbstractArray{T},
                                    soil_lrunoff::AbstractArray{M},
                                    soil_outflux_f::AbstractArray{T},
                                    soil_perc::AbstractArray{M},
                                    soil_freeze_depth::AbstractArray{M},
                                    soil_agtop_cover::AbstractArray{T},
                                    soil_beta_soil::AbstractArray{M},
                                    soil_NO3::AbstractArray{M},
                                    soil_n_leaching::AbstractArray{T},
                                    soil_layer_depth::AbstractArray{T},
                                    kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers, anion_excl, NPERCO = kernel_params

    @unpack soil_infil, soil_infil_litter, percthres = lpjmlparams

    freewater = zero(T)
    soil_srunoff[cell] = zero(T)
    soil_outflux_f[cell] = zero(T)
    soil_n_leaching[cell] = zero(T)
    for l in 1:soil_layers
        freewater += soil_w_fw[l, cell]
        if soil_w[l, cell] + soil_ice_available[l, cell] / soil_whcs[l, cell] > one(T)
            freewater += (soil_w[l, cell] +
                          soil_ice_available[l, cell] / soil_whcs[l, cell] - one(T)) *
                         soil_whcs[l, cell]
        end
        soil_lrunoff[l, cell] = zero(T)
        soil_w_influx[l, cell] = zero(T)
        soil_w_outflux[l, cell] = zero(T)
    end

    soil_infil *= (one(T) + soil_agtop_cover[cell] * soil_infil_litter)
    influx = zero(T)
    iter = 0
    # Iterative slug infiltration + redistribution; hard cap prevents non-convergent loops.
    while (infil[cell] > 1.0f-5 || freewater > 1.0f-5) && iter < 500
        iter += 1
        NO3perc_ly = zero(T)
        freewater = zero(T)
        # Process infiltration in bounded slugs for numerical stability.
        slug = min(4, infil[cell])
        infil[cell] -= slug

        # Calculate influx to first soil layer
        if one(T) - (soil_w[1, cell] * soil_whcs[1, cell] + soil_w_fw[1, cell] + soil_ice_available[1, cell] + soil_ice_free[1, cell]) / (soil_wsats[1, cell] - soil_wpwps[1, cell]) >= 0
            influx = slug * ((1 - (soil_w[1, cell] * soil_whcs[1, cell] + soil_w_fw[1, cell] + soil_ice_available[1, cell] + soil_ice_free[1, cell]) / (soil_wsats[1, cell] - soil_wpwps[1, cell])) ^ (1 / soil_infil))
            soil_w_influx[1, cell] += influx
        else
            influx = zero(T)
            soil_w_influx[1, cell] += influx
        end
        srunoff = slug-influx
        soil_srunoff[cell] += slug - influx # surface runoff used for leaching

        for l in 1:soil_layers
            lrunoff = zero(T)

            # Nitrate percolated from the layer above enters this layer even
            # when water does not continue percolating out of it today.
            soil_NO3[l, cell] += NO3perc_ly
            NO3perc_ly = zero(T)

            soil_w[l, cell] += (soil_w_fw[l, cell] + influx) / soil_whcs[l, cell]
            soil_w_fw[l, cell] = zero(T)

            # Handle lateral runoff of water above saturation
            liquid_capacity = max(
                soil_layer_depth[l] - soil_freeze_depth[l, cell], zero(T),
            ) * (soil_wsat[l, cell] - soil_wpwp[l, cell])
            if (soil_w[l, cell] * soil_whcs[l, cell]) > liquid_capacity
                grunoff = (soil_w[l, cell] * soil_whcs[l, cell]) - liquid_capacity
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end

            # Additional saturation check
            if (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) > soil_wsats[l, cell]
                grunoff = (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) - soil_wsats[l, cell]
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end

            # Percolation from layer l to l+1 (or to outflux at bottom layer).
            if (soil_w[l, cell] - percthres) > (1.0f-5 / soil_whcs[l, cell])
                # Calculate hydraulic conductivity
                ice_fraction = soil_wsats[l, cell] > eps(T) ?
                               clamp(soil_ice[l, cell] / soil_wsats[l, cell], zero(T), one(T)) : zero(T)
                ice_impedance = T(10) ^ (-T(6) * ice_fraction)
                HC = soil_Ks[l, cell] * ((soil_w[l, cell] * soil_whcs[l, cell] + soil_wpwps[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) / soil_wsats[l, cell])^soil_beta_soil[l, cell] * ice_impedance
                # Calculate time constant
                TT = ((soil_w[l, cell] - percthres) * soil_whcs[l, cell]) / HC
                # Calculate percolation amount
                perc = ((soil_w[l, cell] - percthres) * soil_whcs[l, cell]) * (1 - exp(-24 / TT))
                # Correction of percolation for water content of the following layer
                if l < soil_layers
                    saturation_factor = 1 - (soil_w[l+1, cell] * soil_whcs[l+1, cell] + soil_w_fw[l+1, cell] + soil_ice_available[l+1, cell] + soil_ice_free[l+1, cell]) / (soil_wsats[l+1, cell] - soil_wpwps[l+1, cell])
                    if saturation_factor < 0
                        perc = zero(T)
                    else
                        perc *= sqrt(saturation_factor)
                    end
                else
                    saturation_factor = 1 - (soil_w[l, cell] * soil_whcs[l, cell] + soil_w_fw[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) / (soil_wsats[l, cell] - soil_wpwps[l, cell])
                    if saturation_factor < 0
                        perc = zero(T)
                    else
                        perc *= sqrt(saturation_factor)
                    end
                end

                soil_w[l, cell] -= perc / soil_whcs[l, cell]

                if soil_w[l, cell] < 0
                    perc += soil_w[l, cell] * soil_whcs[l, cell]
                    soil_w[l, cell] = zero(T)
                end

                if l == soil_layers
                    soil_outflux_f[cell] += perc
                    soil_w_outflux[l, cell] += perc
                else
                    influx = perc
                    soil_w_influx[l+1, cell] += perc
                    soil_w_outflux[l, cell] += perc
                end

                concNO3_mobile = zero(T)
                # determination of nitrate concentration in mobile water
                w_mobile = perc + srunoff + lrunoff
                if w_mobile > 1.0e-7
                    ww = -w_mobile / ((1 - anion_excl) * soil_wsats[l, cell])
                    vno3 = soil_NO3[l, cell] * (1 - exp(ww))
                    concNO3_mobile = max(vno3 / w_mobile, zero(0))
                end
                # Surface runoff NO3 can be added here if a dedicated pool is introduced.
                srunoff = zero(T)
                NO3lat = zero(T)
                if l == 1
                    NO3lat = NPERCO * concNO3_mobile * lrunoff
                else
                    NO3lat = concNO3_mobile * lrunoff
                end
                NO3lat = min(NO3lat, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3lat
                soil_n_leaching[cell] += NO3lat

                # nitrate percolating from this layer
                NO3perc_ly = concNO3_mobile * perc
                NO3perc_ly = min(NO3perc_ly, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3perc_ly
                if l == soil_layers
                    soil_n_leaching[cell] += NO3perc_ly
                end
            end
        end
    end

    for l in 1:soil_layers
        # Lateral runoff can draw from water already stored in a layer when
        # freezing reduces its liquid pore capacity, so it must be included in
        # the absolute layer-stock update as well as in the outgoing ledger.
        soil_perc[l, cell] = soil_w_influx[l, cell] -
                             soil_w_outflux[l, cell] -
                             soil_lrunoff[l, cell]
    end

end
