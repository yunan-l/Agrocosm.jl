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
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            thermalparams::SoilThermalParams{T} = SoilThermalParams{T}(),
) where {T <: AbstractFloat}
    surface_litter_interception!(soil, prec, crop.fluxes.water.interception)
    transfer_heat = !irrigation && snowmelt !== nothing && air_temperature !== nothing
    if transfer_heat
        infil_perc!(
            soil, prec, snowmelt, air_temperature;
            lpjmlparams = lpjmlparams,
            thermalparams = thermalparams,
        )
    else
        infil_perc!(soil; lpjmlparams = lpjmlparams)
    end

    if !irrigation
        launch_custom!(
            add_layer_flux_kernel!,
            soil.water.storage,
            size(soil.water.storage, 2),
            soil.water.percolation,
            size(soil.water.storage, 1),
        )
    end
    if transfer_heat
        # LPJmL may reconcile temperature every two infiltration iterations.
        # Agrocosm preserves the same water/energy ledger but applies it once
        # after the GPU column kernel, avoiding device synchronization inside
        # the iterative hydrology loop.
        apply_percolation_enthalpy!(soil; thermalparams = thermalparams)
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
        launch_custom!(
            reset_irrigated_storage_kernel!,
            soil.water.storage,
            size(soil.water.storage, 2),
            soil.water.field_capacity,
            soil.properties.layer_depth,
            size(soil.water.storage, 1),
        )
    else
        launch_custom!(
            remove_evapotranspiration_kernel!,
            soil.water.storage,
            size(soil.water.storage, 2),
            crop.fluxes.water.transpiration_layer,
            soil.water.evaporation,
            size(soil.water.storage, 1),
        )
    end
    partition_soil_water_ice!(soil)

    return nothing
end

@kernel inbounds = true function add_layer_flux_kernel!(
    storage::AbstractMatrix{T},
    flux::AbstractMatrix{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        storage[layer, cell] += flux[layer, cell]
    end
end

@kernel inbounds = true function remove_evapotranspiration_kernel!(
    storage::AbstractMatrix{T},
    transpiration::AbstractMatrix{T},
    evaporation::AbstractMatrix{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        storage[layer, cell] -= transpiration[layer, cell] + evaporation[layer, cell]
    end
end

@kernel inbounds = true function reset_irrigated_storage_kernel!(
    storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        storage[layer, cell] = field_capacity[layer, cell] * layer_depth[layer]
    end
end
