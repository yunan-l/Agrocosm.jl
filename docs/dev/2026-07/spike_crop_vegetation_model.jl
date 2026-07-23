# Phase 3 integration spike: a minimal coupled crop model on CPU.
#
# Injects the crop photosynthesis + crop stomatal conductance into Terrarium's
# VegetationCarbon (other slots keep PALADYN defaults), inside a LandModel, and takes a
# coupled timestep — validating the crop physiology runs end-to-end and produces GPP.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_vegetation_model.jl

using Agrocosm
using Terrarium

arch = CPU()
grid = ColumnGrid(arch, ExponentialSpacing(Δz_max = 1.0, N = 20))

# Only the photosynthesis + stomatal-conductance slots are crop here. Swapping additional
# slots (e.g. carbon_dynamics = CropCarbonDynamics(...)) currently fails because Terrarium's
# PALADYN vegetation processes are concretely coupled to each other's types (e.g.
# PALADYNAutotrophicRespiration dispatches on PALADYNCarbonDynamics). Assembling the full crop
# physiology needs a dedicated crop VegetationModel (plan Phase 5), not slot-swapping.
vegetation = VegetationCarbon(
    eltype(grid);
    photosynthesis = CropPhotosynthesis(eltype(grid)),
    stomatal_conductance = CropStomatalConductance(eltype(grid)),
)

hydrology = SoilHydrology(eltype(grid), RichardsEq())
soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
land = LandModel(grid; soil, vegetation)

initializers = (
    temperature = 20.0,
    saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1 * z),
    carbon_vegetation = 0.5,
)
integrator = initialize(land; initializers)

set!(integrator.state.windspeed, 1.0)
set!(integrator.state.specific_humidity, 5.0e-3)
set!(integrator.state.surface_shortwave_down, 400.0)
set!(integrator.state.leaf_area_index, 3.0)

timestep!(integrator, 60.0)

gpp = interior(integrator.state.gross_primary_production)[1, 1, 1]
An = interior(integrator.state.net_assimilation)[1, 1, 1]
Rd = interior(integrator.state.leaf_respiration)[1, 1, 1]
λc = interior(integrator.state.leaf_to_air_co2_ratio)[1, 1, 1]
gc = interior(integrator.state.canopy_water_conductance)[1, 1, 1]

println("SPIKE OK")
println("  photosynthesis        = ", typeof(land.vegetation.photosynthesis).name.name)
println("  stomatal_conductance  = ", typeof(land.vegetation.stomatal_conductance).name.name)
println("  leaf_to_air_co2_ratio = ", λc)
println("  canopy_water_conduct. = ", gc, " m/s")
println("  gross_primary_prod.   = ", gpp, " kgC/m^2/s")
println("  net_assimilation      = ", An, " gC/m^2/s")
println("  leaf_respiration      = ", Rd, " gC/m^2/s")

@assert land.vegetation.photosynthesis isa CropPhotosynthesis
@assert land.vegetation.stomatal_conductance isa CropStomatalConductance
@assert all(isfinite, (gpp, An, Rd, λc, gc))
@assert gpp ≥ 0 "GPP should be non-negative under active conditions"
println("SPIKE ASSERTIONS PASSED")
