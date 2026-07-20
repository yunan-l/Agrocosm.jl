"""Plant carbon stocks and daily carbon flux state."""
mutable struct CropCarbon{A, M}
    biomass::A
    leaf::A
    root::A
    pool::A
    storage::A
    initial_organs::A
    organs::M
    yield::A
    npp::A
    respiration::A
    temperature_response::A
end

function init_crop_carbon(cell_size::Int, device; carbon_pools::Int = 4)
    float_state() = device(zeros(Float32, cell_size))
    initial_organs = Float32[8.0, 0.0113804, 0.0, 11.9886196]

    return CropCarbon(
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        device(initial_organs),
        device(zeros(Float32, carbon_pools, cell_size)),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
    )
end
