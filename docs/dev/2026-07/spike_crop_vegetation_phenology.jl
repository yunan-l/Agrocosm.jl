# Phase 5 spike: a minimal crop VegetationModel driven by prognostic phenological heat units.
#
# Validates that CropVegetation assembles into a Terrarium LandModel and that the heat-unit
# prognostic drives the crop LAI: (1) heat units accumulate over a short integration, and (2) a
# mid-season heat-unit state produces a positive, growing-season LAI and hence positive GPP.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_vegetation_phenology.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# NoFlow soil hydrology (avoids the stiff Richards integration) — the heat-unit → LAI → GPP
# mechanism does not depend on soil water transport for this demonstration.
soil = SoilEnergyWaterCarbon(eltype(grid))
vegetation = CropVegetation(eltype(grid))
land = LandModel(grid; soil, vegetation)

integrator = initialize(land; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)          # warm → heat units accumulate
set!(integrator.state.surface_shortwave_down, 400.0)  # light for photosynthesis

hu0 = interior(integrator.state.phenological_heat_units)[1, 1, 1]
lai0 = interior(integrator.state.leaf_area_index)[1, 1, 1]

# (1) Integrate a few steps: heat units should accumulate above the base temperature.
run!(integrator; steps = 20)
hu1 = interior(integrator.state.phenological_heat_units)[1, 1, 1]

# (2) Jump to a mid-season heat-unit state and recompute auxiliaries; LAI and GPP should be positive.
requirement = land.vegetation.phenology_dynamics.heat_unit_requirement
set!(integrator.state.phenological_heat_units, 0.5 * requirement)   # fphu = 0.5
Terrarium.compute_auxiliary!(integrator.state, integrator.model)

fphu = interior(integrator.state.phenology_heat_unit_fraction)[1, 1, 1]
lai = interior(integrator.state.leaf_area_index)[1, 1, 1]
gpp = interior(integrator.state.gross_primary_production)[1, 1, 1]
β = interior(integrator.state.soil_moisture_limiting_factor)[1, 1, 1]

println("SPIKE OK")
println("  vegetation             = ", typeof(land.vegetation).name.name)
println("  heat units: ", hu0, " -> ", hu1, " (accumulated over 20 steps)")
println("  mid-season fphu        = ", fphu)
println("  leaf_area_index        = ", lai)
println("  soil_moisture_factor β = ", β)
println("  gross_primary_prod.    = ", gpp, " kgC/m^2/s")

@assert land.vegetation isa CropVegetation
@assert hu1 > hu0 "heat units should accumulate at 25 °C"
@assert fphu ≈ 0.5 rtol = 1e-6
@assert lai > 0 "mid-season LAI should be positive"
@assert isfinite(β) "soil-moisture limiting factor should be finite"
@assert gpp > 0 "mid-season GPP should be positive"
println("SPIKE ASSERTIONS PASSED")
