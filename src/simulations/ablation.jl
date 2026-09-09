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
# It does not decide the order. Leaf temperature is solved inside the sub-daily
# loop and the sterility accumulator is filled per sub-step, so both later steps
# rest on the first for physical reasons rather than editorial ones, and the
# cumulative sequence below is the only one in which each rung is a superset of
# the one before it.
#
# The sink is the one step whose prerequisite is weaker than its position:
# it needs the sub-daily loop, not organ temperature. Running it *without*
# organ temperature is not an ordering variant, it is a deliberate off-ladder
# comparison - the accumulator then integrates duration at air temperature
# instead of leaf temperature, which separates the sink mechanism from the
# leaf-air departure that triggers it. `docs/07_ablation_framework.md` explains
# why that cell is load-bearing: on these cells organ temperature's entire yield
# effect arrives through the sink, so the two cannot be attributed apart without
# it.

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

Rung zero is `:daily`, which is not a step but the baseline the ladder departs
from: a daily time step, daily-mean air temperature everywhere, no reproductive
sink, so the harvest index is a constant and yield is set entirely by assimilated
carbon.

It was called `:ggcm` and that name overclaimed. This is *this model with three
processes off*, not a GGCM: it still carries LPJmL's physiology, and it carries
this project's own departures from LPJmL - `senescent_leaf_release = 1`, and the
grain-carbon protection in `carbon_allocation`, which has no switch at all.
Whether rung zero actually reproduces the ensemble's underestimation is an
empirical claim to be checked against GGCMI, not something the ladder
establishes by construction. What the ladder does establish is internal: the
rungs differ in exactly these three fields and nothing else.

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
        # `requires` is the prerequisite the constructor enforces, which for the
        # sink is the sub-daily loop rather than the rung immediately below it.
        # The ladder still places it last because the rungs are cumulative.
        :reproductive_sink, :reproductive_sink, :subdaily,
        "accumulate flowering-window heat exposure and cap the harvest index by grain set",
        "off is bitwise the organ-temperature rung; on with sterility_rate=0 is too",
    ),
)

"""
    ablation_rungs()

Names of the ladder's rungs, baseline first: `(:daily, :subdaily, ...)`.
"""
ablation_rungs() = (:daily, map(step -> step.name, ABLATION_LADDER)...)

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
    # On for every step up to and including `rung`, off above it. `:daily` is
    # below the first step, so nothing is on.
    cut = rung === :daily ? 0 : findfirst(step -> step.name === rung, ABLATION_LADDER)
    settings = (
        step.field => index <= cut
        for (index, step) in enumerate(ABLATION_LADDER)
    )
    return (; settings..., kwargs...)
end

"""
    ablation_air_driven_sink_configuration(; kwargs...)

The off-ladder cell that runs the reproductive sink on **air** temperature: the
sub-daily loop and the sink on, organ temperature off.

Comparing it against the `:reproductive_sink` rung is what separates the sink
mechanism from the leaf-air departure that triggers it. Both runs have the same
sink, the same threshold and the same sub-daily resolution; only the temperature
the accumulator integrates differs. It is also the closest thing in this model
to an air-temperature-driven GGCM sterility function, so the difference is a
statement about the baseline as much as about our own module.

Deliberately not a rung: it is not a structural deficiency of the baseline being
repaired, and folding it into the cumulative ladder would let the headline
decomposition absorb it.
"""
function ablation_air_driven_sink_configuration(; kwargs...)
    owned = map(step -> step.field, ABLATION_LADDER)
    for key in keys(kwargs)
        key in owned && throw(ArgumentError(
            "$key is set by this configuration; it is a fixed comparison cell",
        ))
    end
    return (; subdaily_photosynthesis = true, organ_temperature = false,
            reproductive_sink = true, kwargs...)
end

"""
    ablation_daily_assimilation_sink_configuration(; organ_temperature = true, kwargs...)

The off-ladder cell that runs the reproductive sink on a **daily** assimilation
kernel: sub-daily photosynthesis off, the standalone heat-exposure pass on, the
sink on.

This is the configuration the ladder's own measurements argue for, and the
argument is in `docs/07_ablation_framework.md`. Three findings make it:

  - The exposure accumulator reads only daily state, so it never needed the
    assimilation loop it currently sits inside. Reproducing it externally
    matches the kernel to within 0.04-3.29% at the five gate cells.
  - Sub-daily *assimilation* contributes only -0.1 to -2.3% of the
    daily-invisible event response, while the sink contributes 69-99% of it. The
    expensive half of the sub-daily loop is not the half carrying the signal.
  - Sub-daily assimilation under-assimilates, thinning the canopy and thereby
    overheating the leaf, so the yield deficit against the observational
    reference and the size of the event response are one coupled defect rather
    than two results. The perturbation-B exposure response survives on the daily
    trajectory and at rice is more than twice as large, because a canopy at 39 C
    has already saturated the sterility logistic.

`organ_temperature` is the one switch this cell exposes, because leaf against air
is the comparison that separates the sink mechanism from the departure that
triggers it - the same comparison `ablation_air_driven_sink_configuration` makes
on the sub-daily kernel, and it must stay expressible here too. It defaults to
the leaf, which is the physically right choice and the one the 50-84% exposure
loss at air temperature justifies.

Deliberately not a rung. The ladder is a sequence of structural deficiencies of
the baseline being repaired one at a time; this is a claim about which of two
repairs to keep, and folding it in would let the headline decomposition absorb
the comparison it exists to make.
"""
function ablation_daily_assimilation_sink_configuration(;
    organ_temperature::Bool = true, kwargs...,
)
    owned = (map(step -> step.field, ABLATION_LADDER)..., :subdaily_heat_exposure)
    for key in keys(kwargs)
        key in owned && throw(ArgumentError(
            "$key is set by this configuration; it is a fixed comparison cell " *
            "apart from organ_temperature, which is a named argument",
        ))
    end
    return (; subdaily_photosynthesis = false, subdaily_heat_exposure = true,
            organ_temperature, reproductive_sink = true, kwargs...)
end

"""
    ablation_daily_statistic_sink_configuration(; kwargs...)

The off-ladder floor: the reproductive sink driven by a closed-form duration
above the threshold, computed from the daily mean and the daily range alone -
no sub-step loop, no canopy energy balance, air temperature.

It exists to keep one objection answerable. The paper's control claims a
daily-mean kernel cannot see a perturbation that widens the diurnal range while
holding the mean, and for *this model's* daily kernel that is exactly true. It
is not true of the models the paper positions itself against: they read `tasmax`
and `tasmin`, so widening the range raises their maximum, and a sterility
criterion built on that maximum is not inert. This cell measures what such a
criterion recovers instead of assuming it recovers nothing, which makes the
comparison a result rather than a definition.

Compared upward against `ablation_daily_assimilation_sink_configuration` it
isolates what the canopy energy balance adds over air temperature; compared
against `:daily` it bounds how much of the response a conventional model could
have reached with the forcing it already reads.

The threshold is hard here, where the sub-daily paths smooth it with a 1 C
logistic, and that difference is not incidental - see
`docs/07_ablation_framework.md`, where the smoothing is measured to supply the
majority of accumulated exposure at three of five cells. Reading this cell as
"the closed form loses information" without holding the threshold sharpness
fixed would attribute that smoothing to the integration method.
"""
function ablation_daily_statistic_sink_configuration(; kwargs...)
    owned = (map(step -> step.field, ABLATION_LADDER)...,
             :subdaily_heat_exposure, :daily_statistic_exposure)
    for key in keys(kwargs)
        key in owned && throw(ArgumentError(
            "$key is set by this configuration; it is a fixed comparison cell",
        ))
    end
    return (; subdaily_photosynthesis = false, subdaily_heat_exposure = false,
            daily_statistic_exposure = true, organ_temperature = false,
            reproductive_sink = true, kwargs...)
end

"""
    ablation_terminal_heat_configuration(; organ_temperature = true,
                                         daily_statistic = false, kwargs...)

The off-ladder cell that carries terminal heat and **not** the grain-set sink:
heat damage to grain filling, alone.

Isolating it is the point. The two mechanisms limit different quantities in
sequence, and the harvest index multiplies them, so a run carrying both cannot
say which one a loss came from - and at the hot-wheat cell the answer is not
obvious in advance, because the sink is inert there while the filling window
holds 75.8 exposure hours above 30 C against 0.4 above 35 C. Bounding
`filling_rate` against the observational reference also requires the sink off,
or the bound absorbs grain-set damage at the same cell.

`daily_statistic = true` swaps the sub-daily leaf-temperature pass for the
closed-form daily-statistic path, giving the same floor comparison the sink
cells have. `organ_temperature` is then necessarily off and passing it is an
error rather than a silent override.

`reproductive_sink = true` turns both mechanisms on, which is the complete
reproductive configuration rather than an ablation cell - useful for the
production run and for showing that the two losses compose multiplicatively
rather than one masking the other. It defaults off because the cell's purpose
is isolation.
"""
function ablation_terminal_heat_configuration(; organ_temperature::Bool = true,
                                              daily_statistic::Bool = false,
                                              reproductive_sink::Bool = false,
                                              kwargs...)
    owned = (map(step -> step.field, ABLATION_LADDER)...,
             :subdaily_heat_exposure, :daily_statistic_exposure, :terminal_heat)
    for key in keys(kwargs)
        key in owned && throw(ArgumentError(
            "$key is set by this configuration; it is a fixed comparison cell",
        ))
    end
    daily_statistic && organ_temperature && throw(ArgumentError(
        "the daily-statistic closed form has no sub-steps to solve a canopy " *
        "energy balance at; pass organ_temperature = false with it",
    ))
    return (; subdaily_photosynthesis = false,
            subdaily_heat_exposure = !daily_statistic,
            daily_statistic_exposure = daily_statistic,
            organ_temperature, reproductive_sink, terminal_heat = true,
            kwargs...)
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
over the seasons that COMPLETED. It is reported alongside yield because the two
answer different questions, and a rung that changes one without the other is
telling you which mechanism moved: agronomy puts grain harvest index near
0.4-0.5 for the major cereals and oilseeds, so a plausible yield reached with an
implausible harvest index is right for the wrong reason.

**Both averages divide by seasons completed, not seasons that yielded grain.**
That distinction is the whole point of the metric on an extreme-event ladder,
and getting it wrong makes the table blind to the outcome it exists to measure.
A destroyed season emits zero to `output.crop.yield`
(`src/processes/crop/harvesting.jl`), so dividing by the count of positive
yields removes it from the numerator AND the denominator: two good seasons and
one good plus one destroyed both report the same `yield_per_season`, and a
heatwave that annihilates a crop registers as no change at all. `harvest_date`
is written on the annual emission path whether or not grain was set, so it is
the honest denominator; `harvesting_year` is the flag that separates "ended with
grain" from "ended with nothing".

`finite` is not decoration. A rung that silently produces non-finite carbon is
the failure mode an ablation table would otherwise report as a number.
"""
function ablation_metrics(simulation::CropSimulation)
    crop = simulation.output.crop
    host(field) = Array(getfield(crop, field))
    season_yield = vec(host(:yield))
    season_above = vec(host(:harvest_aboveground_carbon))
    # A season COMPLETED if it reached an emission, destroyed or not. Rows with
    # no season at all have neither a harvest date nor above-ground carbon.
    harvest_date = vec(Array(simulation.output.calendar.harvest_date))
    completed = count(>(0), harvest_date)
    with_grain = count(>(0), season_yield)
    # Destroyed seasons belong in the harvest-index mean at zero, for the same
    # reason they belong in the yield denominator.
    ended = [
        (grain, above)
        for (grain, above) in zip(season_yield, season_above)
        if above > 0
    ]
    gpp = host(:gpp)
    scale = GRAMS_PER_M2_TO_TONNES_PER_HECTARE / CARBON_FRACTION_OF_DRY_MATTER
    total = sum(filter(>(0), season_yield); init = 0.0) * scale
    return (
        seasons = completed,
        seasons_with_grain = with_grain,
        seasons_destroyed = completed - with_grain,
        yield_total = total,
        yield_per_season = completed == 0 ? 0.0 : total / completed,
        harvest_index = isempty(ended) ? 0.0 :
            sum(grain / above for (grain, above) in ended) / length(ended),
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
`:daily`. Both are differences in t ha-1, signed, so a rung that *raises* yield
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
