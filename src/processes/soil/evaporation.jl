"""
evaporation!(pet_eeq, crop, soil)

Compute layer-wise bare-soil evaporation constrained by near-surface water.
"""
function evaporation!(pet_eeq::AbstractArray{T},
                      crop::Crop,
                      soil::Soil;
                      lpjmlparams::LPJmLParams = lpjmlparams

) where {T <: AbstractFloat}

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(evaporation_kernel!,
               pet_eeq,
               crop.canopy.fpar,
               crop.water.transpiration_layer,
               crop.water.canopy_wet,
               soil.water.storage,
               soil.water.wilting_storage,
               soil.water.holding_capacity_storage,
               soil.water.evaporation,
               soil.properties.surface_litter_cover,
               soil.properties.layer_depth,
               kernel_params)

end

@kernel inbounds = true function evaporation_kernel!(
                                     pet_eeq::AbstractArray{T},
                                     crop_fpar::AbstractArray{T},
                                     crop_trans_layer::AbstractArray{M},
                                     crop_canopy_wet::AbstractArray{T},
                                     soil_swc::AbstractArray{M},
                                     soil_wpwps::AbstractArray{M},
                                     soil_whcs::AbstractArray{M},
                                     soil_evap::AbstractArray{M},
                                     soil_agtop_cover::AbstractArray{T},
                                     soil_layer_depth::AbstractArray{T},
                                     kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers = kernel_params

    @unpack PRIESTLEY_TAYLOR = lpjmlparams  # Priestley-Taylor coefficient

    soildepth_evap = lpjmlparams.soildepth_evap

    evap_energy = pet_eeq[cell] * PRIESTLEY_TAYLOR * max(1 - crop_fpar[cell], 0.05)
    # evap_litter = pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - sum(crop_trans_layer[:, cell])

    crop_trans_layer_sum = zero(T)
    for l in 1:soil_layers
        crop_trans_layer_sum += crop_trans_layer[l, cell]
    end

    evap_ratio = zero(T)
    if evap_energy > 1.0f-5 && (pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - crop_trans_layer_sum) > 1.0f-5
        # w_evap is water content in soildepth_evap that can evaporate
        w_evap = zero(T)
        whcs_evap = zero(T)

        for l in 1:soil_layers
            if soildepth_evap > 0
                fraction = min(1, soildepth_evap / soil_layer_depth[l])
                w_evap += (soil_swc[l, cell] - soil_wpwps[l, cell] - crop_trans_layer[l, cell]) * fraction
                whcs_evap += soil_whcs[l, cell] * fraction
                soildepth_evap -= soil_layer_depth[l]
            end
        end

        evap_soil = evap_energy / (1 + exp(5 - 10 * w_evap / whcs_evap)) * max(0.05, (1 - soil_agtop_cover[cell]))
        if w_evap > 0
            evap_ratio = evap_soil / w_evap
        else
            evap_ratio = zero(T)
        end
    end

    soildepth_evap = lpjmlparams.soildepth_evap
    for l in 1:soil_layers
        if soildepth_evap > 0
            fraction = min(1, soildepth_evap / soil_layer_depth[l])
            soil_evap[l, cell] = (soil_swc[l, cell] - soil_wpwps[l, cell]) * evap_ratio * fraction
            soildepth_evap -= soil_layer_depth[l]
        end
    end

end
