"""
Mineral and organic soil nitrogen pools, inputs, and decomposition fluxes.

The shift arrays are fixed normalized post-spin-up distributions; the
`litter_to_*` arrays store the actual daily layer-resolved retained-N fluxes.
"""
mutable struct SoilNitrogen{A, L, M}
    nitrate::M
    ammonium::M
    litter::L
    decomposed_litter::L
    fast::M
    slow::M
    decomposed_fast::M
    decomposed_slow::M
    input::L
    shift_fast::M
    shift_slow::M
    litter_to_fast::M
    litter_to_slow::M
    litter_response::A
    mineralization::M
    immobilization::M
    nitrification::M
    n2o_nitrification::M
    denitrification::M
    n2o_denitrification::M
    n2_denitrification::M
    volatilization::A
    leaching::A
end

init_soil_nitrogen(cell_size::Int, device; kwargs...) =
    init_soil_nitrogen(Float32, cell_size, device; kwargs...)
function init_soil_nitrogen(::Type{T}, cell_size::Int, device;
                            litter_layers::Int = 3,
                            soil_layers::Int = 5) where {T <: AbstractFloat}
    litter_state() = device(zeros(T, litter_layers, cell_size))
    layer_state() = device(zeros(T, soil_layers, cell_size))
    return SoilNitrogen(
        layer_state(), layer_state(),
        litter_state(), litter_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        litter_state(), layer_state(), layer_state(), layer_state(), layer_state(),
        device(zeros(T, litter_layers)),
        layer_state(), layer_state(), layer_state(), layer_state(),
        layer_state(), layer_state(), layer_state(),
        device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
    )
end
