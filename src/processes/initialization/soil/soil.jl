"""Process-grouped soil state with CPU/GPU-compatible array leaves."""
mutable struct Soil{P, W, T, C, N, D, G, L, S}
    properties::P
    water::W
    thermal::T
    carbon::C
    nitrogen::N
    decomposition::D
    management::G
    surface_litter::L
    snow::S
end

init_soil(cell_size::Int, soildepth::AbstractArray{T}, device; kwargs...) where {T <: AbstractFloat} =
    init_soil(T, cell_size, soildepth, device; kwargs...)
function init_soil(::Type{T},
                   cell_size::Int,
                   soildepth::AbstractArray{S},
                   device;
                   litc_layers::Int = 3,
                   soil_layers::Int = 5) where {T <: AbstractFloat, S <: AbstractFloat}
    return Soil(
        init_soil_properties(T, cell_size, soildepth, device),
        init_soil_water(T, cell_size, device; soil_layers = soil_layers),
        init_soil_thermal(T, cell_size, device; soil_layers = soil_layers),
        init_soil_carbon(T, cell_size, device;
                         litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_nitrogen(T, cell_size, device;
                           litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_decomposition(T, cell_size, device; soil_layers = soil_layers),
        init_soil_management(T, cell_size, device; litter_layers = litc_layers),
        init_soil_surface_litter(T, cell_size, device),
        init_soil_snow(T, cell_size, device),
    )
end
