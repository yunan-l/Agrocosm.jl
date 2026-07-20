"""Current-day crop photosynthetic capacities and limiting factors."""
mutable struct CropPhotosynthesisAuxiliary{A}
    potential_vcmax::A
    vcmax::A
    nitrogen_limitation::A
    lambda::A
    temperature_stress::A
end

function init_crop_photosynthesis_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_auxiliary() = device(zeros(T, cell_size))
    return CropPhotosynthesisAuxiliary(ntuple(_ -> float_auxiliary(), 5)...)
end
