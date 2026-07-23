# Phase 2 spike: run a Terrarium coupled `LandModel` column (soil + surface energy
# balance + surface hydrology + atmosphere + vegetation carbon) on CPU.
#
# Purpose: validate that the full coupled land stack Agrocosm will build its crop
# model on computes on this machine, and inspect the exposed state variables
# (soil temperature/moisture, surface fluxes) that the crop processes will couple
# to in Phase 3.
#
# The stack is exercised with a single coupled timestep, matching Terrarium's
# shipped `examples/simulations/land_column.jl`. A sustained multi-step Richards
# integration is stiff at Δt = 60 s under a saturated initial column (the
# saturation prognostic diverges); choosing a stable timestepper/Δt for the crop
# model is a Phase 5 concern. Multi-day *soil* integration is separately
# validated by spike_soil_column_cpu.jl (SoilModel, 3-day run).
#
# Run from the Agrocosm project environment (Terrarium is a dev dependency):
#   julia --project=. docs/dev/2026-07/spike_land_column_cpu.jl

using Terrarium

arch = CPU()

# Soil column with 30 exponentially spaced layers.
grid = ColumnGrid(arch, ExponentialSpacing(Δz_max = 1.0, N = 30))

# Soil hydrology: Richards equation with van Genuchten retention.
swrc = VanGenuchten(α = 2.0, n = 2.0)
hydraulic_properties = ConstantSoilHydraulics(
    eltype(grid); swrc, unsat_hydraulic_cond = UnsatKVanGenuchten(eltype(grid))
)
hydrology = SoilHydrology(eltype(grid), RichardsEq(); hydraulic_properties)

# Coupled soil energy + water + carbon, and a vegetation-carbon component.
soil = SoilEnergyWaterCarbon(eltype(grid); hydrology)
vegetation = VegetationCarbon(eltype(grid))

land = LandModel(grid; soil, vegetation)

# Variably saturated column with a water table near 5 m depth.
initializers = (
    temperature = 15.0,
    saturation_water_ice = (x, z) -> min(1, 0.5 - 0.1 * z),
    carbon_vegetation = 0.5,
)
integrator = initialize(land; initializers)

# Prescribe near-surface atmospheric inputs.
set!(integrator.state.windspeed, 1.0)          # m/s
set!(integrator.state.specific_humidity, 1.0e-4)  # kg/kg

Δt = 60.0
timestep!(integrator, Δt)   # one coupled step (atmosphere→soil→veg→hydrology→SEB)

T = interior(integrator.state.temperature)[1, 1, :]
sat = interior(integrator.state.saturation_water_ice)[1, 1, :]
lhf = interior(integrator.state.latent_heat_flux)[1, 1, 1]
ghf = interior(integrator.state.ground_heat_flux)[1, 1, 1]
gpp = interior(integrator.state.gross_primary_production)[1, 1, 1]
lai = interior(integrator.state.leaf_area_index)[1, 1, 1]

println("SPIKE OK")
println("  n soil layers        = ", length(T))
println("  current_time (s)     = ", current_time(integrator))
println("  soil T[1:3] (°C)     = ", T[1:3])
println("  saturation[1:3]      = ", sat[1:3])
println("  latent_heat_flux     = ", lhf)
println("  ground_heat_flux     = ", ghf)
println("  gross_primary_prod   = ", gpp)
println("  leaf_area_index      = ", lai)
println("  # state variables    = ", length(propertynames(integrator.state)))

for (name, val) in (
        ("soil temperature", T), ("saturation", sat),
        ("latent_heat_flux", [lhf]), ("ground_heat_flux", [ghf]),
        ("gross_primary_production", [gpp]), ("leaf_area_index", [lai]),
    )
    @assert all(isfinite, val) "non-finite $name"
end
println("SPIKE ASSERTIONS PASSED")
