# Phase 5 spike: a crop VegetationModel driven by prognostic phenological heat units and a prognostic
# carbon pool.
#
# Validates that CropVegetation assembles into a Terrarium LandModel and that (1) the heat-unit
# fraction drives a developed growing-season LAI and hence positive GPP and NPP, and (2) the crop
# biomass accumulates NPP over the integration (the carbon loop is closed).
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_vegetation_phenology.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# NoFlow soil hydrology (avoids the stiff Richards integration); the heat-unit → LAI → GPP → biomass
# mechanism does not depend on soil water transport for this demonstration.
soil = SoilEnergyWaterCarbon(eltype(grid))
vegetation = CropVegetation(eltype(grid))
land = LandModel(grid; soil, vegetation)

integrator = initialize(land; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)          # warm → heat units accumulate
set!(integrator.state.surface_shortwave_down, 400.0)  # light for photosynthesis

# Jump to a mid-season heat-unit state so the canopy is developed, then diagnose.
requirement = land.vegetation.phenology_dynamics.heat_unit_requirement
set!(integrator.state.phenological_heat_units, 0.5 * requirement)   # fphu = 0.5
Terrarium.compute_auxiliary!(integrator.state, integrator.model)

fphu = interior(integrator.state.phenology_heat_unit_fraction)[1, 1, 1]
lai = interior(integrator.state.leaf_area_index)[1, 1, 1]
gpp = interior(integrator.state.gross_primary_production)[1, 1, 1]
npp = interior(integrator.state.net_primary_production)[1, 1, 1]
biomass0 = interior(integrator.state.crop_biomass)[1, 1, 1]
hu0 = interior(integrator.state.phenological_heat_units)[1, 1, 1]

# Integrate: crop biomass should accumulate NPP and heat units should keep rising.
run!(integrator; steps = 20)
biomass1 = interior(integrator.state.crop_biomass)[1, 1, 1]
hu1 = interior(integrator.state.phenological_heat_units)[1, 1, 1]

println("SPIKE OK")
println("  vegetation             = ", typeof(land.vegetation).name.name)
println("  mid-season fphu        = ", fphu)
println("  leaf_area_index        = ", lai)
println("  gross_primary_prod.    = ", gpp, " kgC/m^2/s")
println("  net_primary_prod.      = ", npp, " kgC/m^2/s")
println("  crop_biomass: ", biomass0, " -> ", biomass1, " kgC/m^2 (accumulated over 20 steps)")
println("  heat units:   ", hu0, " -> ", hu1)

@assert land.vegetation isa CropVegetation
@assert fphu ≈ 0.5 rtol = 1e-6
@assert lai > 0 "mid-season LAI should be positive"
@assert gpp > 0 "mid-season GPP should be positive"
@assert npp > 0 "mid-season NPP should be positive"
@assert biomass1 > biomass0 "crop biomass should accumulate NPP"
@assert hu1 > hu0 "heat units should keep accumulating at 25 °C"
println("SPIKE ASSERTIONS PASSED")
