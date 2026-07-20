"""Crop-management inputs and grid-cell location state."""
mutable struct ManagedLand{A}
    manure::A
    fertilizer::A
    residue_fraction::A
    latitude::A
end

init_managed_land(cell_size::Int, device) = init_managed_land(Float32, cell_size, device)
function init_managed_land(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return ManagedLand(ntuple(_ -> cell_state(), 4)...)
end
