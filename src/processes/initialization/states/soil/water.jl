"""Layered soil water stocks, hydraulic properties, and daily fluxes."""
mutable struct SoilWater{A, M}
    storage::M
    evaporation::M
    relative_content::M
    free_water::M
    wilting_fraction::M
    wilting_storage::M
    field_capacity::M
    saturation_fraction::M
    saturation_storage::M
    beta::M
    holding_capacity_fraction::M
    holding_capacity_storage::M
    saturated_conductivity::M
    influx::M
    outflux::M
    surface_runoff::A
    lateral_runoff::M
    bottom_drainage::A
    infiltration::A
    percolation::M
end

function init_soil_water(cell_size::Int, device; soil_layers::Int = 5)
    layer_state() = device(zeros(Float32, soil_layers, cell_size))
    cell_state() = device(zeros(Float32, cell_size))
    return SoilWater(
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        cell_state(),
        layer_state(),
        cell_state(),
        cell_state(),
        layer_state(),
    )
end
