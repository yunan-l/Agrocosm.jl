# Phase 5 spike: prognostic soil carbon biogeochemistry in the soil `biogeochem` slot.
#
# Validates that CropSoilBiogeochemistry assembles into a Terrarium SoilModel (feeding density_soc
# into the soil organic fraction) and that the soil carbon pools decompose over the integration with
# positive heterotrophic respiration, driven by the soil temperature/moisture decomposition response.
#
# Run: julia --project=. docs/dev/2026-07/spike_soil_biogeochemistry.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), Float64, ExponentialSpacing(N = 10))

# Soil with the dynamic crop carbon biogeochemistry replacing the constant-density default.
soil = SoilEnergyWaterCarbon(eltype(grid); biogeochem = CropSoilBiogeochemistry(eltype(grid)))
initializer = SoilInitializer(
    eltype(grid);
    energy = QuasiThermalSteadyState(eltype(grid), T₀ = 15.0),   # warm soil → active decomposition
    hydrology = ConstantSaturation(eltype(grid), sat = 0.6),
)
model = SoilModel(grid; soil, initializer)

integrator = initialize(model; boundary_conditions = PrescribedSurfaceTemperature(:T_ub, 15.0))

fast0 = sum(interior(integrator.state.fast_carbon))
slow0 = sum(interior(integrator.state.slow_carbon))
litter0 = sum(interior(integrator.state.litter_carbon))
nh4_0 = sum(interior(integrator.state.soil_ammonium))
no3_0 = sum(interior(integrator.state.soil_nitrate))

run!(integrator; steps = 50)

fast1 = sum(interior(integrator.state.fast_carbon))
slow1 = sum(interior(integrator.state.slow_carbon))
litter1 = sum(interior(integrator.state.litter_carbon))
het = sum(interior(integrator.state.heterotrophic_respiration))
mineralization = sum(interior(integrator.state.net_mineralization))
nh4_1 = sum(interior(integrator.state.soil_ammonium))
no3_1 = sum(interior(integrator.state.soil_nitrate))

println("SPIKE OK")
println("  litter carbon: ", litter0, " -> ", litter1)
println("  fast carbon:   ", fast0, " -> ", fast1)
println("  slow carbon:   ", slow0, " -> ", slow1)
println("  heterotrophic respiration (sum) = ", het, " kgC/m^3/s")
println("  net mineralization (sum)        = ", mineralization, " kgN/m^3/s")
println("  soil ammonium: ", nh4_0, " -> ", nh4_1)
println("  soil nitrate:  ", no3_0, " -> ", no3_1)

@assert soil.biogeochem isa CropSoilBiogeochemistry
@assert all(isfinite, interior(integrator.state.temperature)) "soil temperature stayed finite"
@assert litter1 < litter0 "litter carbon should decompose"
@assert het > 0 "heterotrophic respiration should be positive"
# Total soil carbon decreases (no input yet); the loss balances the respired carbon.
@assert (litter1 + fast1 + slow1) < (litter0 + fast0 + slow0) "total soil carbon should decline without input"
@assert mineralization > 0 "decomposition should mineralize nitrogen"
@assert all(isfinite, interior(integrator.state.soil_nitrate)) "mineral nitrogen stayed finite"
@assert no3_1 > no3_0 "nitrification should build up soil nitrate"
println("SPIKE ASSERTIONS PASSED")
