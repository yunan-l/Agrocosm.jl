"""Static or slowly varying soil physical and chemical properties."""
mutable struct SoilProperties{A, M}
    sand_fraction::M
    clay_fraction::M
    ph::A
    layer_depth::A
end

function init_soil_properties(cell_size::Int, soildepth, device)
    return SoilProperties(
        device(zeros(Float32, 1, cell_size)),
        device(zeros(Float32, 1, cell_size)),
        device(zeros(Float32, cell_size)),
        device(soildepth),
    )
end
