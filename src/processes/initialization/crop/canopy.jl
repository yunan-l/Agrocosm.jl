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

function init_crop_canopy(cell_size::Int, device)
    float_state() = device(zeros(Float32, cell_size))

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
