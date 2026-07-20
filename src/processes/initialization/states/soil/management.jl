const SURFACE_LITTER = 1
const INCORPORATED_LITTER = 2
const ROOT_LITTER = 3

"""Soil management operators and daily internal litter-routing diagnostics."""
mutable struct SoilManagement{M, A}
    tillage_fraction::M
    tillage_carbon::A
    tillage_nitrogen::A
    bioturbation_carbon::A
    bioturbation_nitrogen::A
end

function init_soil_management(cell_size::Int, device; litter_layers::Int = 3)
    cell_state() = device(zeros(Float32, cell_size))
    return SoilManagement(
        device(zeros(Float32, litter_layers, litter_layers)),
        cell_state(), cell_state(), cell_state(), cell_state(),
    )
end
