"""Layered soil temperature, enthalpy, phase state, and thermal properties."""
mutable struct SoilThermal{A, B, M}
    temperature::M
    enthalpy::M
    frozen_fraction::M
    freeze_depth::M
    heat_capacity_frozen::M
    heat_capacity_unfrozen::M
    latent_heat::M
    conductivity_frozen::M
    conductivity_unfrozen::M
    water_reference::M
    surface_energy_flux::A
    energy_residual::A
    initialized::B
    diffusivity_0::A
    diffusivity_15::A
end

function init_soil_thermal(cell_size::Int, device; soil_layers::Int = 5)
    layer_state() = device(zeros(Float32, soil_layers, cell_size))
    cell_state() = device(zeros(Float32, cell_size))
    return SoilThermal(
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
        cell_state(),
        device(zeros(Bool, cell_size)),
        cell_state(),
        cell_state(),
    )
end
