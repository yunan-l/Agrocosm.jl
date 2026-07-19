"""Mineral and organic soil nitrogen pools, inputs, and decomposition fluxes."""
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

function init_soil_nitrogen(cell_size::Int, device;
                            litter_layers::Int = 3,
                            soil_layers::Int = 5)
    litter_state() = device(zeros(Float32, litter_layers, cell_size))
    layer_state() = device(zeros(Float32, soil_layers, cell_size))
    return SoilNitrogen(
        layer_state(), layer_state(),
        litter_state(), litter_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        litter_state(), layer_state(), layer_state(),
        device(zeros(Float32, litter_layers)),
        layer_state(), layer_state(), layer_state(), layer_state(),
        layer_state(), layer_state(), layer_state(),
        device(zeros(Float32, cell_size)),
        device(zeros(Float32, cell_size)),
    )
end
