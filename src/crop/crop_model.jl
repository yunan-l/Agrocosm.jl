# Convenience constructor for a managed-crop land model: a Terrarium `LandModel` assembled from the
# crop vegetation (`CropVegetation`) and a soil with the crop carbon-nitrogen biogeochemistry
# (`CropSoilBiogeochemistry`), configured for a crop functional type from the 12-CFT registry. This
# is the top-level entry point for the crop model — `CropModel(grid, crop_pft(:maize))`.

"""
    $(TYPEDSIGNATURES)

Build a managed-crop `LandModel` on `grid` for the crop functional type `pft` (from the CFT registry,
e.g. `crop_pft(:maize)`): the crop vegetation, and a soil whose biogeochemistry is the crop C–N
scheme. Extra keyword arguments are forwarded to `LandModel` (e.g. `timestepper`, `atmosphere`).
"""
function CropModel(
        grid, pft::PftParameters;
        soil_hydrology = SoilHydrology(eltype(grid)),
        soil_biogeochemistry = CropSoilBiogeochemistry(eltype(grid)),
        vegetation = CropVegetation(eltype(grid), pft),
        kwargs...,
    )
    soil = SoilEnergyWaterCarbon(eltype(grid); hydrology = soil_hydrology, biogeochem = soil_biogeochemistry)
    return LandModel(grid; soil, vegetation, kwargs...)
end

"""
    $(TYPEDSIGNATURES)

Build a managed-crop `LandModel` for the named or numbered crop functional type `crop` (default
temperate cereals). See [`crop_pft`](@ref) for the registry.
"""
CropModel(grid; crop = 1, kwargs...) = CropModel(grid, crop_pft(crop); kwargs...)
