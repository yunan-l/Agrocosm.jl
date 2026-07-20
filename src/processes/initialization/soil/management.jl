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

init_soil_management(cell_size::Int, device; kwargs...) =
    init_soil_management(Float32, cell_size, device; kwargs...)
function init_soil_management(::Type{T}, cell_size::Int, device;
                              litter_layers::Int = 3) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return SoilManagement(
        device(zeros(T, litter_layers, litter_layers)),
        cell_state(), cell_state(), cell_state(), cell_state(),
    )
end
