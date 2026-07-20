"""Crop-management inputs and grid-cell location state."""
mutable struct ManagedLand{A}
    manure::A
    fertilizer::A
    residue_fraction::A
    latitude::A
end

function init_managed_land(cell_size::Int, device)
    cell_state() = device(zeros(Float32, cell_size))
    return ManagedLand(ntuple(_ -> cell_state(), 4)...)
end
