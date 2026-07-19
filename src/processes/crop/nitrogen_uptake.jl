"""
nuptake_crop!(crop, PFT, soil)

Compute root uptake of mineral nitrogen from soil NH4/NO3 pools.
"""
function nuptake_crop!(crop::Crop,
                       PFT::PftParameters,
                       soil::Soil;
                       lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(
        nuptake_crop_kernel!,
        crop.nitrogen,
        crop.nuptake,
        crop.leafn,
        crop.leafc,
        crop.rootn,
        crop.rootc,
        crop.ndemand_leaf,
        crop.ndemand_tot,
        crop.vscal,
        crop.rootdist,
        crop.isgrowing,
        soil.w,
        soil.wsat,
        soil.NO3,
        soil.NH4,
        soil.layer_depth,
        soil.temp,
        PFT,
        kernel_params
    )
  
end

@kernel inbounds = true function nuptake_crop_kernel!(
                                      crop_nitrogen::AbstractArray{T},
                                      crop_nuptake::AbstractArray{T},
                                      crop_leafn::AbstractArray{T},
                                      crop_leafc::AbstractArray{T},
                                      crop_rootn::AbstractArray{T},
                                      crop_rootc::AbstractArray{T},
                                      crop_ndemand_leaf::AbstractArray{T},
                                      crop_ndemand_tot::AbstractArray{T},
                                      crop_vscal::AbstractArray{T},
                                      crop_rootdist::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      soil_w::AbstractArray{M},
                                      soil_wsat::AbstractArray{M},
                                      soil_NO3::AbstractArray{M},
                                      soil_NH4::AbstractArray{M},
                                      soil_layer_depth::AbstractArray{T},
                                      soil_temp::AbstractArray{M},
                                      PFT::PftParameters,
                                      kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)
    
    @unpack lpjmlparams, soil_layers = kernel_params

    @unpack T_0, T_m, T_r = lpjmlparams
    @unpack ncleaf, knstore, no3_uptake, nh4_uptake = PFT

    if crop_isgrowing[cell] == 1
        crop_nuptake[cell] = zero(T)

        mobile_carbon = crop_leafc[cell] + crop_rootc[cell]
        NCplant = mobile_carbon > zero(T) ?
                  (crop_leafn[cell] + crop_rootn[cell]) / mobile_carbon : T(ncleaf.low)
        nc_reference = T(2) / (one(T) / T(ncleaf.low) + one(T) / T(ncleaf.high))
        f_NCplant = clamp(
            (NCplant - T(ncleaf.high)) / (nc_reference - T(ncleaf.high)),
            zero(T),
            one(T),
        )

        leaf_nc = crop_leafc[cell] > zero(T) ?
                  crop_leafn[cell] / crop_leafc[cell] : zero(T)
        total_potential_uptake = zero(T)

        if leaf_nc < T(ncleaf.high) * (one(T) + T(knstore))
            # First pass: independent potential NO3 and NH4 uptake per layer.
            for l in 1:soil_layers
                wscaler = soil_w[l, cell] > T(1e-7) ? one(T) : zero(T)
                temp_response = max(
                    (soil_temp[l, cell] - T(T_0)) *
                    (T(2) * T(T_m) - T(T_0) - soil_temp[l, cell]) /
                    (T(T_r) - T(T_0)) /
                    (T(2) * T(T_m) - T(T_0) - T(T_r)),
                    zero(T),
                )
                root_factor = temp_response * f_NCplant * crop_rootc[cell] *
                              crop_rootdist[l] / T(1000)

                no3_available = max(zero(T), soil_NO3[l, cell])
                if no3_available > zero(T)
                    no3_saturation = no3_available * wscaler /
                                     (no3_available * wscaler + T(no3_uptake.Km) *
                                      soil_wsat[l, cell] * soil_layer_depth[l] / T(1000))
                    no3_potential = T(no3_uptake.vmax) *
                                    (T(no3_uptake.kmin) + no3_saturation) * root_factor
                    total_potential_uptake += min(no3_potential, no3_available)
                end

                nh4_available = max(zero(T), soil_NH4[l, cell])
                if nh4_available > zero(T)
                    nh4_saturation = nh4_available * wscaler /
                                     (nh4_available * wscaler + T(nh4_uptake.Km) *
                                      soil_wsat[l, cell] * soil_layer_depth[l] / T(1000))
                    nh4_potential = T(nh4_uptake.vmax) *
                                    (T(nh4_uptake.kmin) + nh4_saturation) * root_factor
                    total_potential_uptake += min(nh4_potential, nh4_available)
                end
            end
        end

        remaining_demand = max(zero(T), crop_ndemand_tot[cell] - crop_nitrogen[cell])
        n_uptake = min(total_potential_uptake, remaining_demand)

        if n_uptake > zero(T) && total_potential_uptake > zero(T)
            uptake_scale = n_uptake / total_potential_uptake

            # Second pass: remove exactly the accepted uptake from each pool.
            for l in 1:soil_layers
                wscaler = soil_w[l, cell] > T(1e-7) ? one(T) : zero(T)
                temp_response = max(
                    (soil_temp[l, cell] - T(T_0)) *
                    (T(2) * T(T_m) - T(T_0) - soil_temp[l, cell]) /
                    (T(T_r) - T(T_0)) /
                    (T(2) * T(T_m) - T(T_0) - T(T_r)),
                    zero(T),
                )
                root_factor = temp_response * f_NCplant * crop_rootc[cell] *
                              crop_rootdist[l] / T(1000)

                no3_available = max(zero(T), soil_NO3[l, cell])
                if no3_available > zero(T)
                    no3_saturation = no3_available * wscaler /
                                     (no3_available * wscaler + T(no3_uptake.Km) *
                                      soil_wsat[l, cell] * soil_layer_depth[l] / T(1000))
                    no3_potential = min(
                        T(no3_uptake.vmax) * (T(no3_uptake.kmin) + no3_saturation) * root_factor,
                        no3_available,
                    )
                    soil_NO3[l, cell] = max(
                        zero(T), soil_NO3[l, cell] - no3_potential * uptake_scale,
                    )
                end

                nh4_available = max(zero(T), soil_NH4[l, cell])
                if nh4_available > zero(T)
                    nh4_saturation = nh4_available * wscaler /
                                     (nh4_available * wscaler + T(nh4_uptake.Km) *
                                      soil_wsat[l, cell] * soil_layer_depth[l] / T(1000))
                    nh4_potential = min(
                        T(nh4_uptake.vmax) * (T(nh4_uptake.kmin) + nh4_saturation) * root_factor,
                        nh4_available,
                    )
                    soil_NH4[l, cell] = max(
                        zero(T), soil_NH4[l, cell] - nh4_potential * uptake_scale,
                    )
                end
            end

            crop_nitrogen[cell] += n_uptake
            crop_nuptake[cell] = n_uptake
        end

        ndemand_leaf_opt = crop_ndemand_leaf[cell]
        if crop_ndemand_tot[cell] > crop_nitrogen[cell]
            crop_ndemand_leaf[cell] = crop_leafn[cell]
            if ndemand_leaf_opt < 1.0e-7
                crop_vscal[cell] = one(T)
            else
                crop_vscal[cell] = min(one(T), (crop_ndemand_leaf[cell] / (ndemand_leaf_opt / (1 + knstore))))
            end
        else
            crop_vscal[cell] = one(T)
        end

    else
        crop_nitrogen[cell] = zero(T)
        crop_nuptake[cell] = zero(T)
        crop_vscal[cell] = zero(T)
    end
end
