# Phase 4 spike: discrete crop management as Oceananigans callbacks on a Terrarium `Simulation`.
#
# Wraps the full maize `CropModel` in a `Simulation`, registers the sowing and harvest events via
# `add_crop_management!`, and runs a compressed season. Sowing establishes a seeded stand (biomass
# jumps from zero) that grows under the continuous crop dynamics; harvest exports the grain and returns
# the residue to the soil litter (a jump in the soil litter carbon) and clears the stand. The discrete
# events are the sanctioned continuous-time exceptions; everything between them is the continuous
# tendencies.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_management.jl

using Agrocosm
using Terrarium
using Oceananigans: add_callback!, IterationInterval

day = Terrarium.seconds_per_day(Float64)

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
# Compressed phenology (small heat-unit requirement) so the crop completes its cycle within the short
# demonstration season; the crop C–N biogeochemistry soil is needed so harvest residue has a litter
# pool to land in.
vegetation = CropVegetation(eltype(grid); phenology_dynamics = CropPhenologyDynamics(eltype(grid); heat_unit_requirement = 60.0))
model = CropModel(grid, crop_pft("maize"); vegetation = vegetation)
integrator = initialize(model; initializers = (temperature = 20.0,))
set!(integrator.state.air_temperature, 25.0)
set!(integrator.state.surface_shortwave_down, 400.0)

calendar = CropCalendar(Float64; sowing_day = 1, harvest_day = 3, residue_fraction = 0.25)
# Continuous fertilizer applied over days 1–2 (a time-distributed input flux, not a discrete dump).
fertilization = CropFertilization(Float64; application_rate = 1.0e-7, nitrate_fraction = 0.5,
    application_start_day = 1, application_end_day = 2)

# Δt = 600 s matches the stable seasonal spike; larger steps overshoot the surface energy balance.
Δt = 600.0
simulation = Simulation(integrator; Δt = Δt, stop_time = 4 * day)
add_crop_management!(simulation, calendar)
add_crop_fertilization!(simulation, fertilization)

# Record the biomass, soil-litter, and soil mineral-N trajectory every 2 h.
history = NamedTuple{(:day, :biomass, :litter, :mineral_n), NTuple{4, Float64}}[]
function record!(sim)
    s = sim.model.state
    push!(history, (
        day = sim.model.clock.time / day,
        biomass = interior(s.crop_biomass)[1, 1, 1],
        litter = sum(interior(s.litter_carbon)),
        mineral_n = sum(interior(s.soil_ammonium)) + sum(interior(s.soil_nitrate)),
    ))
    return nothing
end
add_callback!(simulation, record!, IterationInterval(12); name = :record)

run!(simulation)

biomass = [h.biomass for h in history]
litter = [h.litter for h in history]
mineral_n = [h.mineral_n for h in history]
peak_biomass = maximum(biomass)
before_harvest = last(filter(h -> h.day < calendar.harvest_day, history))
after_harvest = last(history)
during_fertilization = last(filter(h -> h.day < fertilization.application_end_day, history))

println("SPIKE OK")
println("  day │ crop biomass (kgC/m²) │ soil litter C │ soil mineral N")
for h in history
    println("  ", lpad(round(h.day, digits = 1), 4), " │ ", rpad(round(h.biomass, digits = 5), 12),
        " │ ", rpad(round(h.litter, digits = 4), 9), " │ ", round(h.mineral_n, digits = 5))
end
println("  peak biomass                  = ", round(peak_biomass, digits = 5), " kgC/m²")
println("  biomass just before harvest   = ", round(before_harvest.biomass, digits = 5))
println("  biomass at end (post-harvest) = ", round(after_harvest.biomass, digits = 5))
println("  soil litter: before harvest = ", round(before_harvest.litter, digits = 4),
    " → end = ", round(after_harvest.litter, digits = 4))
println("  soil mineral N: start = ", round(first(mineral_n), digits = 5),
    " → during fertilization = ", round(during_fertilization.mineral_n, digits = 5))

@assert peak_biomass > calendar.seed_carbon "the sown stand should grow beyond the seed"
@assert after_harvest.biomass < before_harvest.biomass "harvest should clear the standing biomass"
@assert after_harvest.litter > before_harvest.litter "harvest should return residue to the soil litter"
@assert during_fertilization.mineral_n > first(mineral_n) "fertilizer should raise the soil mineral nitrogen"
@assert all(isfinite, biomass) && all(isfinite, litter) && all(isfinite, mineral_n) "the managed season stayed finite"
println("SPIKE ASSERTIONS PASSED")
