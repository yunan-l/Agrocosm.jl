# Discrete crop-management events — the sanctioned continuous-time exceptions (Terrarium AGENTS.md).
#
# Sowing and harvest are genuine discrete state transitions in the crop lifecycle that cannot be
# expressed as tendencies: at the sowing date a stand is *established* (a seed carbon/nitrogen pool is
# placed and the phenological clock is reset), and at harvest the crop is *removed* (grain is exported
# off-field, the residue is returned to the soil litter, and the stand is cleared). These are
# implemented as host-side state jumps invoked at specified times through Oceananigans `Callback`s on
# a `Simulation`, per the framework's discrete-event exception clause. Everything that admits a
# continuous formulation — heat-unit accumulation, the carbon/nitrogen dynamics, and fertilizer input
# (see `CropFertilization`) — stays in the continuous tendencies and inputs.

"""
    $(TYPEDEF)

Crop-management calendar: the sowing and harvest dates (days from the simulation time origin) and the
seed/residue parameters governing the discrete stand establishment and removal.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropCalendar{NF}
    "Sowing day (days from the simulation time origin)"
    sowing_day::Int = 1
    "Harvest day (days from the simulation time origin)"
    harvest_day::Int = 300
    "Seed carbon establishing the stand at sowing (kgC/m²; LPJmL seed biomass ≈ 20 gC/m²)"
    seed_carbon::NF = 0.02
    "Seed nitrogen establishing the stand at sowing (kgN/m²; seed C:N ≈ 29)"
    seed_nitrogen::NF = 0.02 / 29
    "Fraction of aboveground residue left on the field at harvest (the rest is exported with the grain)"
    residue_fraction::NF = 0.25
end

CropCalendar(::Type{NF}; kwargs...) where {NF} = CropCalendar{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Simulation time (seconds) of management `day`, measured from the time origin.
"""
management_time(day::Integer, ::Type{NF}) where {NF} = NF(day) * Terrarium.seconds_per_day(NF)

# ---- discrete events ----------------------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Sow the crop: establish the stand by seeding the crop carbon and nitrogen pools and resetting the
phenological heat-unit clock to zero (a fresh season). A discrete state jump — the continuous-time
exception clause applies.
"""
function sow!(integrator, calendar::CropCalendar)
    state = integrator.state
    NF = eltype(state.crop_biomass)
    set!(state.phenological_heat_units, zero(NF))
    set!(state.crop_biomass, NF(calendar.seed_carbon))
    set!(state.crop_nitrogen, NF(calendar.seed_nitrogen))
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Harvest the crop: export the grain (the storage organ plus the non-residue fraction of the leaf
carbon/nitrogen) off-field, return the residue (the residue fraction of the leaf plus all of the root
carbon/nitrogen) to the soil litter distributed over the root zone (mass-conserving, since the root
fraction sums to unity), and clear the stand (biomass, nitrogen, and heat units back to zero). Returns
the total exported grain carbon (kgC/m², domain sum) as a harvest-yield diagnostic. A discrete state
jump — the continuous-time exception clause applies.
"""
function harvest!(integrator, calendar::CropCalendar)
    state = integrator.state
    grid = get_grid(integrator.model)
    NF = eltype(state.crop_biomass)
    residue_fraction = NF(calendar.residue_fraction)

    # Yield diagnostic (read before the pools are cleared): grain = storage organ + exported leaf.
    grain_carbon = sum(interior(state.storage_carbon)) +
        (one(NF) - residue_fraction) * sum(interior(state.leaf_carbon))

    # Return the residue to the soil litter / mineral nitrogen over the root zone.
    out = (litter_carbon = state.litter_carbon, soil_ammonium = state.soil_ammonium)
    fields = (
        leaf_carbon = state.leaf_carbon, root_carbon = state.root_carbon,
        leaf_nitrogen = state.leaf_nitrogen, root_nitrogen = state.root_nitrogen,
        root_fraction = state.root_fraction,
    )
    launch!(grid, XYZ, harvest_residue_kernel!, out, fields, residue_fraction)

    # Clear the stand.
    set!(state.crop_biomass, zero(NF))
    set!(state.crop_nitrogen, zero(NF))
    set!(state.phenological_heat_units, zero(NF))
    return grain_carbon
end

@kernel inbounds = true function harvest_residue_kernel!(out, grid, fields, residue_fraction)
    i, j, k = @index(Global, NTuple)
    field_grid = get_field_grid(grid)
    # Distribute the per-area residue over the root zone as a per-volume increment (÷ layer thickness);
    # the aboveground residue (leaf·residue_fraction) is root-weighted alongside the root residue — an
    # approximation that concentrates it near the surface where the roots are densest.
    per_volume = fields.root_fraction[i, j, k] / Δzᵃᵃᶜ(i, j, k, field_grid)
    residue_carbon = (fields.leaf_carbon[i, j] * residue_fraction + fields.root_carbon[i, j]) * per_volume
    residue_nitrogen = (fields.leaf_nitrogen[i, j] * residue_fraction + fields.root_nitrogen[i, j]) * per_volume
    out.litter_carbon[i, j, k] += residue_carbon
    out.soil_ammonium[i, j, k] += residue_nitrogen
end

# ---- Oceananigans callback glue -----------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

Register the crop `calendar`'s sowing and harvest as discrete Oceananigans `Callback`s on the
`simulation` (whose model is a Terrarium `ModelIntegrator`), each scheduled at its management day via
`SpecifiedTimes`. This is the sanctioned mechanism for the crop lifecycle's discrete events; the
continuous crop and soil dynamics run in the timestepper as usual between them.
"""
function add_crop_management!(simulation, calendar::CropCalendar{NF}) where {NF}
    add_callback!(simulation, sim -> sow!(sim.model, calendar),
        SpecifiedTimes(management_time(calendar.sowing_day, NF)); name = :crop_sowing)
    add_callback!(simulation, sim -> harvest!(sim.model, calendar),
        SpecifiedTimes(management_time(calendar.harvest_day, NF)); name = :crop_harvest)
    return simulation
end

# ---- continuous fertilizer application ----------------------------------------------------

"""
    $(TYPEDEF)

Fertilizer application specification: a mineral-nitrogen application rate split between the ammonium
and nitrate pools, applied over a window of the season. Unlike sowing and harvest this is *not* a
discrete event — it is realized as a continuous input flux to the soil biogeochemistry (see
[`fertilize!`](@ref)), so the nitrogen enters through the soil mineral-N tendency and is integrated in
time (a time-distributed input, per the continuous-where-feasible management rule).

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropFertilization{NF}
    "Total mineral-nitrogen application rate while active (kgN/m²/s)"
    application_rate::NF = 0.0
    "Fraction of the application supplied as nitrate (the rest as ammonium)"
    nitrate_fraction::NF = 0.5
    "First day of the application window (days from the simulation time origin)"
    application_start_day::Int = 1
    "Last day of the application window (days from the simulation time origin)"
    application_end_day::Int = 300
end

CropFertilization(::Type{NF}; kwargs...) where {NF} = CropFertilization{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Set the soil biogeochemistry's fertilizer input fluxes from the application `fert` and the current
simulation time: the split application rate while the clock is within the window, zero otherwise. Meant
to be called each step (e.g. via [`add_crop_fertilization!`](@ref)) so the continuous flux tracks the
clock; the soil biogeochemistry integrates it into the mineral-nitrogen pools.
"""
function fertilize!(integrator, fert::CropFertilization)
    state = integrator.state
    NF = eltype(state.fertilizer_ammonium_flux)
    time = state.clock.time
    active = management_time(fert.application_start_day, NF) <= time <
        management_time(fert.application_end_day, NF)
    rate = ifelse(active, NF(fert.application_rate), zero(NF))
    set!(state.fertilizer_ammonium_flux, rate * (one(NF) - NF(fert.nitrate_fraction)))
    set!(state.fertilizer_nitrate_flux, rate * NF(fert.nitrate_fraction))
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Register the continuous fertilizer application `fert` as a per-step Oceananigans `Callback` on the
`simulation`, keeping the soil biogeochemistry's fertilizer input fluxes current with the clock so the
application is integrated continuously over its window.
"""
function add_crop_fertilization!(simulation, fert::CropFertilization)
    add_callback!(simulation, sim -> fertilize!(sim.model, fert),
        IterationInterval(1); name = :crop_fertilization)
    return simulation
end
