"""
    soil_infiltration!(soil, crop, precipitation)

Apply throughfall infiltration and percolation before the daily plant water-stress calculation.
Rainfed and irrigated simulations share this hydrologic and thermal update.
Explicit surface irrigation enters with precipitation before this call;
legacy reset irrigation restores storage later in `soil_evapotranspiration!`.
"""
function soil_infiltration!(soil,
                            crop,
                            prec::AbstractArray{T};
                            snowmelt::Union{Nothing, AbstractArray{T}} = nothing,
                            air_temperature::Union{Nothing, AbstractArray{T}} = nothing,
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            thermalparams::SoilThermalParams{T} = SoilThermalParams{T}(),
) where {T <: AbstractFloat}
    surface_litter_interception!(soil, prec, crop_fluxes(crop).water.interception)
    # The same infiltration, percolation, and enthalpy transport applies to
    # rainfall and explicitly supplied surface irrigation water.
    transfer_heat = snowmelt !== nothing && air_temperature !== nothing
    if transfer_heat
        infil_perc!(
            soil, prec, snowmelt, air_temperature;
            lpjmlparams = lpjmlparams,
            thermalparams = thermalparams,
        )
    else
        infil_perc!(soil; lpjmlparams = lpjmlparams)
    end

    launch_custom!(
        add_layer_flux_kernel!,
        soil_water_prognostic(soil).storage,
        size(soil_water_prognostic(soil).storage, 2),
        soil_water_fluxes(soil).percolation,
        size(soil_water_prognostic(soil).storage, 1),
    )
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
    surface_irrigation!(applied, water_input, state, precipitation, canopy_interception, cft)

Apply enough surface water to refill the upper 500 mm to field capacity when
root-weighted plant-available water falls below 0.8 (C3), 0.7 (C4), or
1.0 (rice). Same-day
liquid throughfall and snowmelt reduce the requested amount. This is a field-level,
unlimited-supply approximation, without conveyance losses or irrigation-
technology-specific physics. Water subsequently follows ordinary infiltration.
"""
function surface_irrigation!(applied::AbstractVector{T},
                             water_input::AbstractVector{T},
                             state::ModelState, precipitation::AbstractVector{T},
                             canopy_interception::AbstractVector{T},
                             cft::CFTParameters) where {T <: AbstractFloat}
    threshold = cft.name == 2 ? one(T) : cft.path == 2 ? T(0.7) : T(0.8)
    launch_1D!(
        surface_irrigation_kernel!, applied, water_input,
        crop_prognostic(state).phenology.is_growing,
        soil_water_auxiliary(state).relative_content,
        crop_root_input(state).distribution,
        soil_water_prognostic(state).storage,
        soil_water_prognostic(state).ice_storage,
        soil_water_auxiliary(state).field_capacity,
        soil_properties(state).layer_depth,
        precipitation, canopy_interception, threshold,
    )
    return nothing
end

@inline function surface_irrigation_amount(
    cell::Integer,
    is_growing::AbstractVector{S},
    relative_content::AbstractMatrix{T},
    root_distribution::AbstractVector{T},
    storage::AbstractMatrix{T},
    ice_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    precipitation::AbstractVector{T},
    canopy_interception::AbstractVector{T},
    threshold::T,
) where {T <: AbstractFloat, S <: Integer}
    requested = zero(T)
    if is_growing[cell] == 1
        root_water = zero(T)
        for layer in 1:5
            root_water += relative_content[layer, cell] * root_distribution[layer]
        end
        if root_water < threshold
            remaining_depth = T(500)
            for layer in 1:5
                depth = layer_depth[layer]
                fraction = min(max(remaining_depth, zero(T)), depth) / depth
                target = field_capacity[layer, cell] * depth
                requested += max(target - storage[layer, cell] - ice_storage[layer, cell], zero(T)) * fraction
                remaining_depth -= depth
            end
            throughfall = max(precipitation[cell] - canopy_interception[cell], zero(T))
            requested = max(requested - throughfall, zero(T))
        end
    end
    return requested
end

@kernel inbounds = true function surface_irrigation_kernel!(
    applied::AbstractVector{T},
    water_input::AbstractVector{T},
    is_growing::AbstractVector{S},
    relative_content::AbstractMatrix{T},
    root_distribution::AbstractVector{T},
    storage::AbstractMatrix{T},
    ice_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    precipitation::AbstractVector{T},
    canopy_interception::AbstractVector{T},
    threshold::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    requested = surface_irrigation_amount(
        cell, is_growing, relative_content, root_distribution,
        storage, ice_storage, field_capacity, layer_depth,
        precipitation, canopy_interception, threshold,
    )
    applied[cell] = requested
    water_input[cell] = precipitation[cell] + requested
end


"""
    soil_evapotranspiration!(soil, crop; irrigation=false)

Remove the current day's layer-resolved transpiration and soil evaporation after
plant water demand and supply have been calculated.
"""
function soil_evapotranspiration!(soil,
                                  crop;
                                  irrigation = false)
    if irrigation
        launch_custom!(
            reset_irrigated_storage_kernel!,
            soil_water_prognostic(soil).storage,
            size(soil_water_prognostic(soil).storage, 2),
            soil_water_prognostic(soil).ice_storage,
            soil_water_auxiliary(soil).field_capacity,
            soil_properties(soil).layer_depth,
            size(soil_water_prognostic(soil).storage, 1),
        )
    else
        launch_custom!(
            remove_evapotranspiration_kernel!,
            soil_water_prognostic(soil).storage,
            size(soil_water_prognostic(soil).storage, 2),
            crop_fluxes(crop).water.transpiration_layer,
            soil_water_fluxes(soil).evaporation,
            size(soil_water_prognostic(soil).storage, 1),
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
    ice_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        target_storage = field_capacity[layer, cell] * layer_depth[layer]
        storage[layer, cell] = max(
            target_storage - ice_storage[layer, cell],
            zero(T),
        )
    end
end
