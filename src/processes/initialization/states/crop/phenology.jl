"""Phenology and crop-calendar state stored as backend arrays."""
mutable struct CropPhenology{A, B, I}
    phu::A
    vdsum::A
    husum::A
    fphu::A
    senescence::B
    senescence_previous::B
    harvesting::B
    harvesting_previous::B
    growing_days::I
    winter_type::B
    is_growing::I
end

function init_crop_phenology(cell_size::Int, device)
    float_state() = device(zeros(Float32, cell_size))
    bool_state(value = false) = device(fill(value, cell_size))

    return CropPhenology(
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        bool_state(),
        bool_state(),
        bool_state(true),
        bool_state(true),
        device(zeros(Int32, cell_size)),
        bool_state(),
        device(zeros(Int32, cell_size)),
    )
end
