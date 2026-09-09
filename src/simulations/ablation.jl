# Controlled within-model ablation.
#
# The attribution argument this model exists to support is not "our model is
# better". It is "we isolated each structural deficiency and measured what it
# contributes", which only works if the comparison is internal: one model, one
# parameter set, processes switched on one at a time. That makes the switch
# inventory a scientific object, not a convenience, and it is declared here so
# that the ladder a paper reports and the ladder the code can actually run are
# the same list.
#
# Two things this file deliberately does NOT do.
#
# It does not re-implement the switches. Every entry names a field of
# `SimulationConfiguration`, and `ablation_configuration` returns keywords for
# the ordinary constructor, so an ablation run and a production run go through
# identical validation. A registry that could describe a configuration the
# constructor rejects would be worse than no registry.
#
# It does not decide the order. The order is forced by the prerequisites, which
# are physical rather than editorial: leaf temperature is solved inside the
# sub-daily loop, and the sterility accumulator reads leaf temperature per
# sub-step. That is why the ladder has exactly one valid sequence, and why
# "which order did you switch them on" is not a degree of freedom a reviewer
# needs to worry about here.

"""
    AblationStep

One switchable process in the controlled ablation ladder.

`field` is the `SimulationConfiguration` field the step sets, `requires` names
the step that must already be on (`nothing` for the first), and `retreat`
records what switching it off is guaranteed to reproduce. `retreat` is
documentation of a verified claim, not an assertion evaluated here; the tests
that establish each one are named in `docs/03_verification_record.md`.
"""
struct AblationStep
    name::Symbol
    field::Symbol
    requires::Union{Nothing, Symbol}
    summary::String
    retreat::String
end

"""
    ABLATION_LADDER

The controlled ablation sequence, in the only order its prerequisites permit.

Rung zero is `:ggcm`, which is not a step but the baseline the ladder departs
from: a daily time step driven by daily-mean air temperature, i.e. what the
GGCM ensemble does. Each subsequent entry switches on one process and leaves
every earlier one on.

`docs/05_reproductive_sink_design.md` and `docs/01_subdaily_design.md` carry the
physics; the point of the table is that the *paper's* ladder and the code's are
literally the same object. It is meant to be reproduced in the supplementary
material, which is where the model description now lives.
"""
const ABLATION_LADDER = (
    AblationStep(
        :subdaily, :subdaily_photosynthesis, nothing,
        "integrate assimilation over the day instead of at the daily mean",
        "off is bitwise the daily kernel; on with steps=1 and :flat is too",
    ),
    AblationStep(
        :organ_temperature, :organ_temperature, :subdaily,
        "solve canopy energy balance per sub-step and drive leaf processes with leaf temperature",
        "off is bitwise the sub-daily rung",
    ),
    AblationStep(
        :reproductive_sink, :reproductive_sink, :organ_temperature,
        "accumulate flowering-window heat exposure and cap the harvest index by grain set",
        "off is bitwise the organ-temperature rung; on with sterility_rate=0 is too",
    ),
)

"""
    ablation_rungs()

Names of the ladder's rungs, baseline first: `(:ggcm, :subdaily, ...)`.
"""
ablation_rungs() = (:ggcm, map(step -> step.name, ABLATION_LADDER)...)

"""
    ablation_step(name) -> AblationStep

The ladder entry called `name`, or an error naming the valid rungs.
"""
function ablation_step(name::Symbol)
    for step in ABLATION_LADDER
        step.name === name && return step
    end
    throw(ArgumentError(
        "unknown ablation step $name; the ladder is $(ablation_rungs())",
    ))
end

"""
    ablation_configuration(rung; kwargs...)

Keyword settings that put a run at ladder rung `rung`, to be splatted into
`initialize_simulation`.

Every process at or below `rung` is on and everything above it is off, so the
returned settings are a complete statement of the run's structure rather than a
patch on whatever the defaults happen to be. That matters for the paper's
central control: Fig 1 claims identical parameters across rungs, and it can
only claim that if the rungs differ in exactly these fields and nothing else.

`kwargs` are passed through unchanged, which is how a rung is combined with the
settings the ladder says nothing about (`subdaily_steps`, `diurnal_shape`,
irrigation, fertiliser). Passing a field the ladder owns is an error rather than
a silent override.

```julia
initialize_simulation(cft, prepared; days,
                      ablation_configuration(:organ_temperature)...)
```
"""
function ablation_configuration(rung::Symbol; kwargs...)
    rungs = ablation_rungs()
    rung in rungs || throw(ArgumentError(
        "unknown ablation rung $rung; the ladder is $rungs",
    ))
    owned = map(step -> step.field, ABLATION_LADDER)
    for key in keys(kwargs)
        key in owned && throw(ArgumentError(
            "$key is set by the ablation rung; pass a different rung instead of " *
            "overriding it, or the run is no longer a point on the ladder",
        ))
    end
    # On for every step up to and including `rung`, off above it. `:ggcm` is
    # below the first step, so nothing is on.
    cut = rung === :ggcm ? 0 : findfirst(step -> step.name === rung, ABLATION_LADDER)
    settings = (
        step.field => index <= cut
        for (index, step) in enumerate(ABLATION_LADDER)
    )
    return (; settings..., kwargs...)
end

"""
    ablation_ladder_settings(; kwargs...)

The whole ladder as `rung => settings` pairs, in order, for a driver that runs
every rung of one case. `kwargs` are shared by every rung, which is what keeps
the comparison controlled.
"""
ablation_ladder_settings(; kwargs...) =
    map(rung -> rung => ablation_configuration(rung; kwargs...), ablation_rungs())

# ---------------------------------------------------------------------------
# What a rung reports
# ---------------------------------------------------------------------------
#
# The ladder is only useful if every rung is summarised the same way, so the
# summary lives here rather than in whichever script happens to drive a run.
# The harvest index in particular has to be defined once: it is the quantity
# that separates "the crop grew less" from "the crop set less grain", which is
# the distinction the reproductive sink exists to make and the one a
# photosynthesis-only model cannot express.

"""Carbon fraction of crop dry matter, the divisor that turns gC into g DM."""
const CARBON_FRACTION_OF_DRY_MATTER = 0.45

"""g m-2 to t ha-1."""
const GRAMS_PER_M2_TO_TONNES_PER_HECTARE = 0.01

"""
    ablation_metrics(simulation) -> NamedTuple

Summarise a finished run for one rung of the ablation ladder.

Aggregated over cells and harvested seasons, which suits the small diagnostic
domains the ladder is exercised on; a global run should reduce
`simulation.output` directly rather than collapse it here.

`harvest_index` is grain carbon over above-ground carbon at harvest, averaged
over harvested seasons. It is reported alongside yield because the two answer
different questions, and a rung that changes one without the other is telling
you which mechanism moved: agronomy puts grain harvest index near 0.4-0.5 for
the major cereals and oilseeds, so a plausible yield reached with an
implausible harvest index is right for the wrong reason.

`finite` is not decoration. A rung that silently produces non-finite carbon is
the failure mode an ablation table would otherwise report as a number.
"""
function ablation_metrics(simulation::CropSimulation)
    crop = simulation.output.crop
    host(field) = Array(getfield(crop, field))
    season_yield = vec(host(:yield))
    season_above = vec(host(:harvest_aboveground_carbon))
    harvested = filter(>(0), season_yield)
    filled = [
        (grain, above)
        for (grain, above) in zip(season_yield, season_above)
        if grain > 0 && above > 0
    ]
    gpp = host(:gpp)
    scale = GRAMS_PER_M2_TO_TONNES_PER_HECTARE / CARBON_FRACTION_OF_DRY_MATTER
    total = isempty(harvested) ? 0.0 : sum(harvested) * scale
    return (
        seasons = length(harvested),
        yield_total = total,
        yield_per_season = isempty(harvested) ? 0.0 : total / length(harvested),
        harvest_index = isempty(filled) ? 0.0 :
            sum(grain / above for (grain, above) in filled) / length(filled),
        season_gpp = sum(gpp),
        peak_lai = maximum(host(:lai); init = 0.0),
        peak_biomass = maximum(host(:biomass); init = 0.0),
        finite = all(isfinite, gpp) && all(isfinite, host(:npp)),
    )
end

"""
    ablation_report(run_rung; shared...)

Run every rung of the ladder and return `(; rows, marginal, total)`.

`run_rung(rung, settings)` builds, runs and returns a finished `CropSimulation`
for one rung; the caller owns input assembly, which is what keeps data loading
out of the model. `shared` keywords go to every rung unchanged - that sharing is
the experiment's control, not a convenience.

`rows` is one `ablation_metrics` per rung with its name attached. `marginal`
gives each rung's change in `yield_per_season` against the rung below it, which
is the decomposition the paper reports; `total` gives its change against
`:ggcm`. Both are differences in t ha-1, signed, so a rung that *raises* yield
is visible rather than folded into a magnitude.
"""
function ablation_report(run_rung; kwargs...)
    rows = [
        (; rung, ablation_metrics(run_rung(rung, settings))...)
        for (rung, settings) in ablation_ladder_settings(; kwargs...)
    ]
    baseline = first(rows).yield_per_season
    marginal = [
        rows[index].rung => rows[index].yield_per_season -
                            rows[index - 1].yield_per_season
        for index in 2:length(rows)
    ]
    total = [
        row.rung => row.yield_per_season - baseline for row in rows[2:end]
    ]
    return (; rows, marginal, total)
end
