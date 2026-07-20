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
    percolation_energy::M
    surface_energy_flux::A
    energy_residual::A
    untracked_water_energy_flux::A
    rain_energy_input::A
    snowmelt_energy_input::A
    lateral_runoff_energy_output::A
    bottom_drainage_energy_output::A
    percolation_energy_residual::A
    initialized::B
    diffusivity_0::A
    diffusivity_15::A
end

init_soil_thermal(cell_size::Int, device; kwargs...) =
    init_soil_thermal(Float32, cell_size, device; kwargs...)
function init_soil_thermal(::Type{T}, cell_size::Int, device;
                           soil_layers::Int = 5) where {T <: AbstractFloat}
    layer_state() = device(zeros(T, soil_layers, cell_size))
    cell_state() = device(zeros(T, cell_size))
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
        layer_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        device(zeros(Bool, cell_size)),
        cell_state(),
        cell_state(),
    )
end
