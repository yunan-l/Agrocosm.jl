"""
    soil_infiltration!(soil, crop, precipitation; irrigation=false)

Apply throughfall infiltration and percolation before the daily plant water-stress calculation.
For rainfed simulations, the resulting layer balance is immediately
added to absolute soil water storage.
"""
function soil_infiltration!(soil::Soil,
                            crop::Crop,
                            prec::AbstractArray{T};
                            irrigation = false
) where {T <: AbstractFloat}
    soil.infil .= prec - crop.intercep
    infil_perc!(soil)

    if !irrigation
        soil.swc .= soil.swc .+ soil.perc
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
        soil.swc .= soil.wfc .* soil.layer_depth
    else
        soil.swc .= soil.swc .- crop.trans_layer .- soil.evap
    end

    return nothing
end
