"""Persistent seasonal crop-water memory."""
mutable struct CropWaterState{A}
    demand_sum::A
    supply_sum::A
    waterlogging_days::A
end

"""Current-day crop-water fluxes."""
mutable struct CropWaterFluxes{A, M}
    transpiration::A
    interception::A
    transpiration_layer::M
end

function init_crop_water_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropWaterState(ntuple(_ -> float_state(), 3)...)
end

function init_crop_water_fluxes(::Type{T}, cell_size::Int, device;
                                soil_layers::Int = 5) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropWaterFluxes(
        float_flux(),
        float_flux(),
        device(zeros(T, soil_layers, cell_size)),
    )
end
