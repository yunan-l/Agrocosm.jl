"""Static or slowly varying soil physical and chemical properties."""
mutable struct SoilProperties{A, M}
    sand_fraction::M
    clay_fraction::M
    ph::A
    layer_depth::A
end

init_soil_properties(cell_size::Int, soildepth, device) =
    init_soil_properties(Float32, cell_size, soildepth, device)
function init_soil_properties(::Type{T}, cell_size::Int, soildepth, device) where {T <: AbstractFloat}
    return SoilProperties(
        device(zeros(T, 1, cell_size)),
        device(zeros(T, 1, cell_size)),
        device(zeros(T, cell_size)),
        device(T.(soildepth)),
    )
end
