"""Environmental response multipliers shared by soil C and N decomposition."""
mutable struct SoilDecomposition{L, M}
    response::M
    litter_response::L
end

function init_soil_decomposition(cell_size::Int, device; soil_layers::Int = 5)
    return SoilDecomposition(
        device(zeros(Float32, soil_layers, cell_size)),
        device(zeros(Float32, 1, cell_size)),
    )
end
