"""Canopy geometry and radiation state stored as backend arrays."""
mutable struct CropCanopy{A}
    lai::A
    flaimax::A
    laimax_adjusted::A
    lai_npp_deficit::A
    phenology_fraction::A
    albedo::A
    fpar::A
    apar::A
end

init_crop_canopy(cell_size::Int, device) = init_crop_canopy(Float32, cell_size, device)
function init_crop_canopy(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))

    return CropCanopy(
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
    )
end
