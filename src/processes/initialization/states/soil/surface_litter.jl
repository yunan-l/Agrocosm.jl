"""Hydrological and thermal state of the above-ground litter layer."""
mutable struct SoilSurfaceLitter{A}
    dry_matter::A
    depth::A
    cover::A
    water_capacity::A
    water_storage::A
    interception::A
    evaporation::A
    temperature::A
    conductivity::A
end

function init_soil_surface_litter(cell_size::Int, device)
    state() = device(zeros(Float32, cell_size))
    return SoilSurfaceLitter(ntuple(_ -> state(), 9)...)
end
