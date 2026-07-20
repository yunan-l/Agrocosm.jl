"""Canopy variables whose previous-day values affect the next daily state."""
mutable struct CropCanopyState{A}
    lai::A
    laimax_adjusted::A
    lai_npp_deficit::A
end

"""Current-day canopy geometry, radiation, and conductance variables."""
mutable struct CropCanopyAuxiliary{A}
    flaimax::A
    phenology_fraction::A
    albedo::A
    fpar::A
    apar::A
    canopy_conductance::A
    canopy_wet::A
end

function init_crop_canopy_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropCanopyState(ntuple(_ -> float_state(), 3)...)
end

function init_crop_canopy_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_auxiliary() = device(zeros(T, cell_size))
    return CropCanopyAuxiliary(
        float_auxiliary(),
        float_auxiliary(),
        float_auxiliary(),
        float_auxiliary(),
        float_auxiliary(),
        float_auxiliary(),
        float_auxiliary(),
    )
end
