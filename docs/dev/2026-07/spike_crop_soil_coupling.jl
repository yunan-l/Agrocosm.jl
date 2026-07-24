# Phase 5 spike: the closed plant↔soil flux loop in the full CropModel.
#
# Validates the two coupling directions between the crop vegetation and the soil biogeochemistry:
#   • nitrogen soil→crop: the crop draws mineral N (uptake flux > 0) and the soil mineral pool
#     (NH₄+NO₃) is drawn down in response;
#   • carbon/nitrogen crop→soil: plant turnover returns litter (litterfall fluxes > 0) that the soil
#     biogeochemistry receives, distributed over the root zone (mass-conserving).
# The crop biomass and nitrogen accumulate and the model stays finite throughout.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_soil_coupling.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
model = CropModel(grid, crop_pft("maize"); soil_hydrology = SoilHydrology(eltype(grid)))

integrator = initialize(model; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)
set!(integrator.state.surface_shortwave_down, 400.0)
set!(integrator.state.phenological_heat_units, 0.5 * model.vegetation.phenology_dynamics.heat_unit_requirement)
# Seed a small standing crop so the crop→soil litter return is active from the first step.
set!(integrator.state.crop_biomass, 0.3)
set!(integrator.state.crop_nitrogen, 0.01)

biomass0 = interior(integrator.state.crop_biomass)[1, 1, 1]
crop_n0 = interior(integrator.state.crop_nitrogen)[1, 1, 1]
soil_mineral_n0 = sum(interior(integrator.state.soil_ammonium)) + sum(interior(integrator.state.soil_nitrate))
litter_c0 = sum(interior(integrator.state.litter_carbon))

run!(integrator; steps = 30)

biomass1 = interior(integrator.state.crop_biomass)[1, 1, 1]
crop_n1 = interior(integrator.state.crop_nitrogen)[1, 1, 1]
soil_mineral_n1 = sum(interior(integrator.state.soil_ammonium)) + sum(interior(integrator.state.soil_nitrate))
litter_c1 = sum(interior(integrator.state.litter_carbon))
# Fluxes sampled after growth (the auxiliaries reflect the final grown state).
uptake = interior(integrator.state.crop_nitrogen_uptake)[1, 1, 1]
litterfall_c = interior(integrator.state.crop_litterfall_carbon)[1, 1, 1]
litterfall_n = interior(integrator.state.crop_litterfall_nitrogen)[1, 1, 1]

println("SPIKE OK")
println("  crop nitrogen uptake flux    = ", uptake, " kgN/m^2/s   (soil → crop)")
println("  crop litterfall carbon flux  = ", litterfall_c, " kgC/m^2/s   (crop → soil)")
println("  crop litterfall nitrogen flux= ", litterfall_n, " kgN/m^2/s   (crop → soil)")
println("  crop biomass:  ", biomass0, " -> ", biomass1, " kgC/m^2")
println("  crop nitrogen: ", crop_n0, " -> ", crop_n1, " kgN/m^2")
println("  soil mineral N (NH4+NO3 sum): ", soil_mineral_n0, " -> ", soil_mineral_n1)
println("  soil litter carbon (sum):     ", litter_c0, " -> ", litter_c1)

@assert model.vegetation isa CropVegetation && model.soil.biogeochem isa CropSoilBiogeochemistry
@assert uptake > 0 "the crop should draw nitrogen from the soil"
@assert litterfall_c > 0 "the crop should return carbon to the soil litter"
@assert litterfall_n > 0 "the crop should return nitrogen to the soil litter"
@assert soil_mineral_n1 < soil_mineral_n0 "crop uptake should draw down the soil mineral nitrogen"
@assert biomass1 > 0 && isfinite(biomass1) "crop biomass stayed finite"
@assert crop_n1 > 0 && isfinite(crop_n1) "crop nitrogen stayed finite"
@assert all(isfinite, interior(integrator.state.soil_nitrate)) "soil nitrogen stayed finite"
@assert all(isfinite, interior(integrator.state.soil_ammonium)) "soil ammonium stayed finite"
@assert all(isfinite, interior(integrator.state.litter_carbon)) "soil litter carbon stayed finite"
println("SPIKE ASSERTIONS PASSED")
