"""Environmental response multipliers shared by soil C and N decomposition."""
mutable struct SoilDecomposition{A, M}
    response::M
    litter_response::M
    layer_scratch_1::M
    layer_scratch_2::M
    surface_scratch_1::A
    surface_scratch_2::A
end

init_soil_decomposition(cell_size::Int, device; kwargs...) =
    init_soil_decomposition(Float32, cell_size, device; kwargs...)
function init_soil_decomposition(::Type{T}, cell_size::Int, device;
                                 soil_layers::Int = 5,
                                 litter_layers::Int = 3) where {T <: AbstractFloat}
    layer_state() = device(zeros(T, soil_layers, cell_size))
    return SoilDecomposition(
        layer_state(),
        device(zeros(T, litter_layers, cell_size)),
        layer_state(),
        layer_state(),
        device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
    )
end
