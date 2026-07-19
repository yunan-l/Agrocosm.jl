"""
    soil_infiltration!(soil, crop, precipitation; irrigation=false)

Apply throughfall infiltration and percolation before the daily plant water-stress calculation.
For rainfed simulations, the resulting layer balance is immediately
added to absolute soil water storage.
"""
function soil_infiltration!(soil::Soil,
                            crop::Crop,
                            prec::AbstractArray{T};
                            irrigation = false,
                            snowmelt::Union{Nothing, AbstractArray{T}} = nothing,
                            air_temperature::Union{Nothing, AbstractArray{T}} = nothing,
) where {T <: AbstractFloat}
    soil.water.infiltration .= prec - crop.water.interception
    surface_litter_interception!(soil)
    transfer_heat = !irrigation && snowmelt !== nothing && air_temperature !== nothing
    if transfer_heat
        infil_perc!(soil, prec, snowmelt, air_temperature)
    else
        infil_perc!(soil)
    end

    if !irrigation
        soil.water.storage .= soil.water.storage .+ soil.water.percolation
    end
    if transfer_heat
        # LPJmL may reconcile temperature every two infiltration iterations.
        # Agrocosm preserves the same water/energy ledger but applies it once
        # after the GPU column kernel, avoiding device synchronization inside
        # the iterative hydrology loop.
        apply_percolation_enthalpy!(soil)
    else
        partition_soil_water_ice!(soil)
    end

    return nothing
end


"""
    soil_evapotranspiration!(soil, crop; irrigation=false)

Remove the current day's layer-resolved transpiration and soil evaporation after
plant water demand and supply have been calculated.
"""
function soil_evapotranspiration!(soil::Soil,
                                  crop::Crop;
                                  irrigation = false)
    if irrigation
        # Preserve the existing idealized irrigation behaviour.
        soil.water.storage .= soil.water.field_capacity .* soil.properties.layer_depth
    else
        soil.water.storage .= soil.water.storage .- crop.water.transpiration_layer .- soil.water.evaporation
    end
    partition_soil_water_ice!(soil)

    return nothing
end
