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

function init_soil(cell_size::Int,
                   soildepth::AbstractArray{T},
                   device;
                   litc_layers::Int = 3,
                   soil_layers::Int = 5) where {T <: AbstractFloat}
    return Soil(
        init_soil_properties(cell_size, soildepth, device),
        init_soil_water(cell_size, device; soil_layers = soil_layers),
        init_soil_thermal(cell_size, device; soil_layers = soil_layers),
        init_soil_carbon(cell_size, device;
                         litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_nitrogen(cell_size, device;
                           litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_decomposition(cell_size, device; soil_layers = soil_layers),
        init_soil_management(device; litter_layers = litc_layers),
        init_soil_surface_litter(cell_size, device),
        init_soil_snow(cell_size, device),
    )
end
