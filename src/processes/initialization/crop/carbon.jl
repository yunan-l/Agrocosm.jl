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

init_crop_carbon(cell_size::Int, device; kwargs...) =
    init_crop_carbon(Float32, cell_size, device; kwargs...)
function init_crop_carbon(::Type{T}, cell_size::Int, device;
                          carbon_pools::Int = 4) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    initial_organs = T[8.0, 0.0113804, 0.0, 11.9886196]

    return CropCarbon(
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        device(initial_organs),
        device(zeros(T, carbon_pools, cell_size)),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
    )
end
