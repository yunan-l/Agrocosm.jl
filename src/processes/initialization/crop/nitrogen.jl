"""Plant nitrogen stocks, demands, management buffers, and stress state."""
mutable struct CropNitrogen{A}
    total::A
    uptake::A
    auto_fertilizer::A
    leaf::A
    root::A
    pool::A
    storage::A
    demand_total::A
    demand_leaf::A
    pending_manure::A
    pending_fertilizer::A
    seed_input::A
    prescribed_manure_input::A
    prescribed_fertilizer_input::A
    harvest_export::A
    stress_sum::A
    stress::A
    deficit::A
end

function init_crop_nitrogen(cell_size::Int, device)
    float_state() = device(zeros(Float32, cell_size))
    return CropNitrogen(ntuple(_ -> float_state(), 18)...)
end
