"""
crop_nitrogen!(crop, PFT, soil, photos_vcmax, temp; auto_fertilizer=true)

Allocate acquired crop nitrogen among leaf, root, storage, and pool compartments.
"""
function crop_nitrogen!(crop::Crop,
                        PFT::PftParameters,
                        soil::Soil,
                        photos_vcmax::AbstractArray{T},
                        temp::AbstractArray{T};
                        auto_fertilizer::Bool = true,
                        lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}

    ndemand_crop!(crop, PFT, photos_vcmax, temp; lpjmlparams = lpjmlparams)
    nuptake_crop!(
        crop, PFT, soil;
        auto_fertilizer = auto_fertilizer,
        lpjmlparams = lpjmlparams,
    )

    allocate_crop_nitrogen!(crop, PFT)

end

"""
    allocate_crop_nitrogen!(crop, PFT)

Redistribute the complete plant nitrogen stock among crop organs. Organ pools
are derived stocks, not daily uptake fluxes, so repeated calls with unchanged
carbon and total nitrogen are idempotent.
"""
function allocate_crop_nitrogen!(crop::Crop,
                                 PFT::PftParameters)

    launch_1D!(crop_nitrogen_kernel!,
               crop.state.nitrogen.total,
               crop.state.phenology.is_growing,
               crop.state.carbon.leaf,
               crop.state.carbon.root,
               crop.state.carbon.storage,
               crop.state.carbon.pool,
               crop.state.nitrogen.leaf,
               crop.state.nitrogen.root,
               crop.state.nitrogen.storage,
               crop.state.nitrogen.pool,
               PFT)

end


@kernel inbounds = true function crop_nitrogen_kernel!(
                                       crop_nitrogen::AbstractArray{T},
                                       crop_isgrowing::AbstractArray{S},
                                       crop_leafc::AbstractArray{T},
                                       crop_rootc::AbstractArray{T},
                                       crop_stoc::AbstractArray{T},
                                       crop_poolc::AbstractArray{T},
                                       crop_leafn::AbstractArray{T},
                                       crop_rootn::AbstractArray{T},
                                       crop_ston::AbstractArray{T},
                                       crop_pooln::AbstractArray{T},
                                       PFT::PftParameters
) where {T <: AbstractFloat, S <: Integer}

     cell = @index(Global)

     @unpack ratio = PFT

     if (crop_isgrowing[cell] == 1) && (crop_nitrogen[cell] > zero(T)) && (crop_leafc[cell] > T(1e-7))
          # LPJmL clears all four organ pools before calling solve(). With those
          # inputs fixed at zero, solve() reduces exactly to these positive
          # carbon-to-target-C:N weights. This form avoids reading stale organ N,
          # uses fewer operations in a GPU kernel, and conserves total plant N.
          leaf_weight = crop_leafc[cell]
          root_weight = crop_rootc[cell] / T(ratio.root)
          sto_weight = crop_stoc[cell] / T(ratio.sto)
          pool_weight = crop_poolc[cell] / T(ratio.pool)
          total_weight = leaf_weight + root_weight + sto_weight + pool_weight

          if total_weight > zero(T)
               scale = crop_nitrogen[cell] / total_weight
               crop_leafn[cell] = leaf_weight * scale
               crop_rootn[cell] = root_weight * scale
               crop_ston[cell] = sto_weight * scale
               crop_pooln[cell] = pool_weight * scale
          else
               crop_leafn[cell] = zero(T)
               crop_rootn[cell] = zero(T)
               crop_ston[cell] = zero(T)
               crop_pooln[cell] = zero(T)
          end
     else
          crop_leafn[cell] = zero(T)
          crop_rootn[cell] = zero(T)
          crop_ston[cell] = zero(T)
          crop_pooln[cell] = zero(T)
     end
end
