"""
Soil carbon inputs, pools, decomposition fluxes, and respiration.

`shift_fast` and `shift_slow` are fixed post-spin-up vertical distributions;
each sums to one per cell. `litter_to_fast` and `litter_to_slow` are the daily
layer-resolved fluxes after applying fastfrac and the retained carbon fraction.
"""
mutable struct SoilCarbon{A, L, M}
    input::L
    litter::L
    decomposed_litter::L
    fast::M
    slow::M
    decomposed_fast::M
    decomposed_slow::M
    shift_fast::M
    shift_slow::M
    litter_to_fast::M
    litter_to_slow::M
    litter_response::A
    heterotrophic_respiration::A
end

function init_soil_carbon(cell_size::Int, device;
                          litter_layers::Int = 3,
                          soil_layers::Int = 5)
    litter_state() = device(zeros(Float32, litter_layers, cell_size))
    layer_state() = device(zeros(Float32, soil_layers, cell_size))
    return SoilCarbon(
        litter_state(), litter_state(), litter_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        device(zeros(Float32, litter_layers)),
        device(zeros(Float32, cell_size)),
    )
end
