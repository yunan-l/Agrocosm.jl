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
               crop.auxiliary.canopy.fpar,
               crop.fluxes.water.transpiration_layer,
               crop.auxiliary.canopy.canopy_wet,
               soil.water.relative_content,
               soil.water.free_water,
               soil.water.holding_capacity_storage,
               soil.water.evaporation,
               soil.surface_litter.cover,
               soil.surface_litter.water_capacity,
               soil.surface_litter.water_storage,
               soil.surface_litter.evaporation,
               soil.properties.layer_depth,
               kernel_params)

end

@kernel inbounds = true function evaporation_kernel!(
                                     pet_eeq::AbstractArray{T},
                                     crop_fpar::AbstractArray{T},
                                     crop_trans_layer::AbstractArray{M},
                                     crop_canopy_wet::AbstractArray{T},
                                     soil_w::AbstractArray{M},
                                     soil_w_fw::AbstractArray{M},
                                     soil_whcs::AbstractArray{M},
                                     soil_evap::AbstractArray{M},
                                     soil_agtop_cover::AbstractArray{T},
                                     litter_water_capacity::AbstractArray{T},
                                     litter_water_storage::AbstractArray{T},
                                     litter_evaporation::AbstractArray{T},
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

    available_evaporation = max(
        pet_eeq[cell] * PRIESTLEY_TAYLOR *
        (one(T) - crop_canopy_wet[cell]) - crop_trans_layer_sum,
        zero(T),
    )
    if litter_water_capacity[cell] > eps(T) && available_evaporation > zero(T)
        litter_wetness = clamp(
            litter_water_storage[cell] / litter_water_capacity[cell],
            zero(T), one(T),
        )
        litter_evaporation[cell] = min(
            evap_energy * litter_wetness * litter_wetness * soil_agtop_cover[cell],
            litter_water_storage[cell],
            available_evaporation,
        )
        litter_water_storage[cell] -= litter_evaporation[cell]
    else
        litter_evaporation[cell] = zero(T)
    end

    evap_ratio = zero(T)
    if evap_energy > 1.0f-5 && (pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - crop_trans_layer_sum) > 1.0f-5
        # w_evap is water content in soildepth_evap that can evaporate
        w_evap = zero(T)
        whcs_evap = zero(T)

        for l in 1:soil_layers
            if soildepth_evap > 0
                fraction = min(1, soildepth_evap / soil_layer_depth[l])
                liquid_above_pwp = soil_w[l, cell] * soil_whcs[l, cell] +
                                   soil_w_fw[l, cell]
                w_evap += max(liquid_above_pwp - crop_trans_layer[l, cell], zero(T)) * fraction
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
            liquid_above_pwp = soil_w[l, cell] * soil_whcs[l, cell] +
                               soil_w_fw[l, cell]
            soil_evap[l, cell] = liquid_above_pwp * evap_ratio * fraction
            soildepth_evap -= soil_layer_depth[l]
        end
    end

end
