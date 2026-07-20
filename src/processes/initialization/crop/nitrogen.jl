"""Persistent plant nitrogen stocks and seasonal memory."""
mutable struct CropNitrogenState{A}
    total::A
    leaf::A
    root::A
    pool::A
    storage::A
    pending_manure::A
    pending_fertilizer::A
    stress_sum::A
end

"""Current-day plant and management nitrogen fluxes."""
mutable struct CropNitrogenFluxes{A}
    uptake::A
    auto_fertilizer::A
    seed_input::A
    prescribed_manure_input::A
    prescribed_fertilizer_input::A
    harvest_export::A
end

function init_crop_nitrogen_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropNitrogenState(ntuple(_ -> float_state(), 8)...)
end

function init_crop_nitrogen_fluxes(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropNitrogenFluxes(ntuple(_ -> float_flux(), 6)...)
end
