"""Soil management operators used by cultivation and residue incorporation."""
mutable struct SoilManagement{M}
    tillage_fraction::M
end

function init_soil_management(device; litter_layers::Int = 3)
    return SoilManagement(device(zeros(Float32, litter_layers, litter_layers)))
end
