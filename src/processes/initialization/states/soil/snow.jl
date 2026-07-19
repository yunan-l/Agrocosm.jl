"""Snow storage, phase-change fluxes, and surface snow diagnostics."""
mutable struct SoilSnow{A}
    pack::A
    melt::A
    sublimation::A
    runoff::A
    height::A
    fraction::A
end

function init_soil_snow(cell_size::Int, device)
    cell_state() = device(zeros(Float32, cell_size))
    return SoilSnow(ntuple(_ -> cell_state(), 6)...)
end
