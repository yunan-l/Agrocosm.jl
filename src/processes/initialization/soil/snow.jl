"""Snow storage, phase-change fluxes, and surface snow diagnostics."""
mutable struct SoilSnow{A}
    pack::A
    melt::A
    sublimation::A
    runoff::A
    height::A
    fraction::A
end

init_soil_snow(cell_size::Int, device) = init_soil_snow(Float32, cell_size, device)
function init_soil_snow(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return SoilSnow(ntuple(_ -> cell_state(), 6)...)
end
