"""Daily crop photosynthesis fluxes, capacities, and stress state."""
mutable struct CropPhotosynthesis{A}
    gross_assimilation::A
    net_assimilation::A
    water_limited_assimilation::A
    leaf_respiration::A
    potential_vmax::A
    vmax::A
    nitrogen_limitation::A
    lambda::A
    temperature_stress::A
end

init_crop_photosynthesis(cell_size::Int, device) =
    init_crop_photosynthesis(Float32, cell_size, device)
function init_crop_photosynthesis(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropPhotosynthesis(ntuple(_ -> float_state(), 9)...)
end
