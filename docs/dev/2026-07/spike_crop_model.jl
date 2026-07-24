# Phase 5 spike: the top-level managed-crop model (CropModel) with both the crop vegetation and the
# soil carbon-nitrogen biogeochemistry active in one Terrarium LandModel.
#
# Validates that CropModel assembles for a CFT and runs end-to-end: the crop produces GPP (from the
# heat-unit-driven canopy) while the soil biogeochemistry cycles carbon and nitrogen.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_model.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# Maize (CFT 3, C4). NoFlow soil hydrology keeps the demonstration off the stiff Richards path.
model = CropModel(grid, crop_pft("maize"); soil_hydrology = SoilHydrology(eltype(grid)))

integrator = initialize(model; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)
set!(integrator.state.surface_shortwave_down, 400.0)
# Develop the canopy (mid-season heat units).
set!(integrator.state.phenological_heat_units, 0.5 * model.vegetation.phenology_dynamics.heat_unit_requirement)

het0 = sum(interior(integrator.state.heterotrophic_respiration))
no3_0 = sum(interior(integrator.state.soil_nitrate))

timestep!(integrator, 60.0)   # one coupled step of the full managed-crop model

gpp = interior(integrator.state.gross_primary_production)[1, 1, 1]
lai = interior(integrator.state.leaf_area_index)[1, 1, 1]
biomass = interior(integrator.state.crop_biomass)[1, 1, 1]
het = sum(interior(integrator.state.heterotrophic_respiration))
no3 = sum(interior(integrator.state.soil_nitrate))
nh4 = sum(interior(integrator.state.soil_ammonium))

println("SPIKE OK")
println("  vegetation / pathway   = C4 maize (", typeof(model.vegetation.photosynthesis).name.name, ")")
println("  leaf_area_index        = ", lai)
println("  gross_primary_prod.    = ", gpp, " kgC/m^2/s")
println("  crop_biomass           = ", biomass, " kgC/m^2")
println("  heterotrophic resp.    = ", het, " kgC/m^3/s")
println("  soil ammonium (sum)    = ", nh4)
println("  soil nitrate (sum)     = ", no3)

@assert model.vegetation isa CropVegetation
@assert model.soil.biogeochem isa CropSoilBiogeochemistry
@assert lai > 0 && gpp > 0 "the crop canopy should assimilate carbon"
@assert het > 0 "the soil should respire carbon"
@assert all(isfinite, interior(integrator.state.temperature)) "soil temperature stayed finite"
@assert all(isfinite, interior(integrator.state.soil_nitrate)) "soil nitrogen stayed finite"
println("SPIKE ASSERTIONS PASSED")
