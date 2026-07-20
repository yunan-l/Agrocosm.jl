"""
    update_surface_litter_properties!(soil; thermalparams=soil_thermal_params)

Update LPJmL-style above-ground litter dry matter, cover, depth, and water
capacity from the first litter-carbon pool. Water exceeding a reduced capacity
is returned to the first soil layer to conserve total water.
"""
function update_surface_litter_properties!(
    soil::Soil;
    thermalparams::SoilThermalParams = soil_thermal_params,
)
    launch_1D!(
        update_surface_litter_properties_kernel!,
        soil.surface_litter.dry_matter,
        soil.carbon.litter,
        soil.surface_litter.depth,
        soil.surface_litter.cover,
        soil.surface_litter.water_capacity,
        soil.surface_litter.water_storage,
        soil.water.storage,
        thermalparams,
    )
    return nothing
end

@kernel inbounds = true function update_surface_litter_properties_kernel!(
    dry_matter::AbstractArray{T},
    litter_carbon::AbstractArray{T},
    depth::AbstractArray{T},
    cover::AbstractArray{T},
    water_capacity::AbstractArray{T},
    water_storage::AbstractArray{T},
    soil_water_storage::AbstractArray{T},
    thermalparams::SoilThermalParams{T},
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack litter_carbon_fraction, litter_bulk_density = thermalparams

    dm = max(litter_carbon[1, cell], zero(T)) / litter_carbon_fraction
    dry_matter[cell] = dm
    depth[cell] = dm / T(1000) / litter_bulk_density
    cover[cell] = one(T) - exp(-T(6e-3) * dm)
    capacity = T(2e-3) * dm
    water_capacity[cell] = capacity

    if water_storage[cell] > capacity
        soil_water_storage[1, cell] += water_storage[cell] - capacity
        water_storage[cell] = capacity
    end
end

"""Fill the surface-litter water store from canopy throughfall."""
function surface_litter_interception!(soil::Soil)
    launch_1D!(
        surface_litter_interception_kernel!,
        soil.water.infiltration,
        soil.surface_litter.cover,
        soil.surface_litter.water_capacity,
        soil.surface_litter.water_storage,
        soil.surface_litter.interception,
    )
    return nothing
end

"""Compute canopy throughfall and surface-litter interception in one cell kernel."""
function surface_litter_interception!(soil::Soil,
                                      precipitation::AbstractVector{T},
                                      canopy_interception::AbstractVector{T}) where {T <: AbstractFloat}
    launch_1D!(
        surface_litter_throughfall_kernel!,
        soil.water.infiltration,
        precipitation,
        canopy_interception,
        soil.surface_litter.cover,
        soil.surface_litter.water_capacity,
        soil.surface_litter.water_storage,
        soil.surface_litter.interception,
    )
    return nothing
end

@kernel inbounds = true function surface_litter_throughfall_kernel!(
    throughfall::AbstractVector{T},
    precipitation::AbstractVector{T},
    canopy_interception::AbstractVector{T},
    cover::AbstractVector{T},
    water_capacity::AbstractVector{T},
    water_storage::AbstractVector{T},
    interception::AbstractVector{T},
) where {T <: AbstractFloat}
    cell = @index(Global)
    incoming = precipitation[cell] - canopy_interception[cell]
    available_capacity = max(water_capacity[cell] - water_storage[cell], zero(T))
    captured = min(available_capacity, max(incoming, zero(T)) * cover[cell])
    water_storage[cell] += captured
    throughfall[cell] = incoming - captured
    interception[cell] = captured
end

@kernel inbounds = true function surface_litter_interception_kernel!(
    throughfall::AbstractArray{T},
    cover::AbstractArray{T},
    water_capacity::AbstractArray{T},
    water_storage::AbstractArray{T},
    interception::AbstractArray{T},
) where {T <: AbstractFloat}
    cell = @index(Global)
    available_capacity = max(water_capacity[cell] - water_storage[cell], zero(T))
    captured = min(available_capacity, max(throughfall[cell], zero(T)) * cover[cell])
    water_storage[cell] += captured
    throughfall[cell] -= captured
    interception[cell] = captured
end
