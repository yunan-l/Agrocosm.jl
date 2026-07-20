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

init_soil_surface_litter(cell_size::Int, device) =
    init_soil_surface_litter(Float32, cell_size, device)
function init_soil_surface_litter(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    state() = device(zeros(T, cell_size))
    return SoilSurfaceLitter(ntuple(_ -> state(), 9)...)
end
