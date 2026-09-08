# Conditional weather-to-yield attribution

This opt-in implementation is isolated from the production checkout. It uses
the pinned source revision's default CFT and process parameters, nitrogen
limitation, fixed prescribed management, existing irrigation and 24-substep
soil heat. It does not fit model parameters or optimize weather.

## Numerical contract

`weather_forcing(climate)` packs prepared model-unit matrices into
`(day, cell, variable)` order: temperature (°C), precipitation (mm/day),
shortwave and net longwave (W/m²), and wind speed (m/s). CO₂ and nitrogen
deposition remain fixed. Calendar indexing is 365-day, including across years.

`enzyme_weather_harvest_gradient` accepts one established crop, its post-sowing
state, consecutive growth days, and the following harvest day. Reverse AD
propagates the state cotangent backward through recomputed blocks. C3 and C4
have separate photosynthesis/lambda dispatch. As in the existing C3 AD
adapter, a Newton-refined lambda is used for derivatives; every case must
match an independent ordinary-production harvest before its gradient is used.

Normal harvest transfers storage carbon before that day's growth. The
fixed-event terminal observable is therefore pre-harvest storage carbon
divided by 0.45 and multiplied by 0.01, yielding t dry matter/ha. This conversion
does not apply to a quantity already reported in dry matter. Harvest-day
weather has zero derivative in this conditional objective: a possible change
in the harvest date is deliberately outside it.

The production `readclimate!` inactive Enzyme rule is retained for existing
parameter objectives. Only the explicit weather-control API activates the five
weather fields. Its forcing-copy kernel supports CPU/CUDA. Full-model weather
AD currently uses CPU `Array` and a single cell; GPU forcing compatibility is
not a claim of full-model GPU AD support.

## What the result does and does not establish

- Gradients are **local, model-conditional sensitivities**, not observational
  proof that weather caused a real-world harvest failure.
- The initial crop/soil state, sowing event, parameters, management, CO₂ and N
  deposition are held fixed. Earlier weather effects on establishment and
  antecedent soil conditions are not included.
- Discrete harvest-date shifts and crop failure are checked with the ordinary
  lifecycle. They must not be interpreted as a smooth fixed-event derivative.
  An unharvested counterfactual is recorded as missing/NaN, not zero yield.
- Reference runs substitute all five weather variables from an explicitly
  selected, calendar-matched reference year (and its preceding year), starting
  after factual sowing. This retains within-day covariance, but does not imply
  that the reference year is an unbiased causal counterfactual.
- `sum(gradient * weather_difference)` is a first-order estimate, not an exact
  or integrated-gradient decomposition. The runner reports its residual
  against the ordinary finite counterfactual, along with event changes.
- Carbon, water, nitrogen and phenology outputs are daily **process
  diagnostics**. Paired weather-window experiments now connect their responses
  to an ordinary finite yield change; they are not additive, causally identified
  process contributions. Formal process-path decomposition is a separate step.
- Attribution cannot recover omitted damage mechanisms (for example a missing
  direct reproductive heat-injury pathway). Only daily mean temperature is
  controlled here, not unobserved hourly peaks or humidity.

## Weather–window–process–yield evidence chain

In addition to full-reference weather runs, the case runner restores reference
weather in **disjoint, fixed calendar windows**, by default 14 days from the first
post-sowing growth day. The last window may be shorter. All five weather channels
are replaced together; every other weather day, the initial state, parameters,
CO₂, deposition and management inputs stay factual. The ordinary model evolves
all states freely and follows any subsequent harvest/failure change. There is no
process override, unlimited-N comparison or parameter optimization.

Window selection is specified before observing the window responses, rather
than testing only days with the largest gradient. For each reference year, all
windows are reported, including zero/negative yield responses. The NetCDF records
the exact replacement mask, 365-day dates, factual PHU fraction, local weather
gradients and ordinary outcomes. A weather difference relative to one reference
year is not by itself proof of a climatological extreme. Reference-year choice
and splice-boundary weather discontinuities remain sensitivity considerations.

The process chain includes:

- temperature response, APAR, canopy conductance and potential/realized Vcmax;
- soil water, transpiration, evaporation and the existing daily water-sufficiency
  multiplier; `rootzone_available_water` is the model's **top-three-layer,
  root-weighted, pre-extraction diagnostic**, not total end-of-day soil water;
- whole-column nitrate/ammonium, mineralization, leaching, volatilization,
  crop N uptake, leaf N and the retained N-supported Vcmax fraction;
- actual daily prescribed fertilizer/manure inputs: the amounts and decision
  rule remain fixed, but weather-driven PHU changes can move the model's second
  application (`fphu > 0.25`). Record this endogenous event shift rather than
  mislabelling every N response as altered soil supply;
- GPP, plant respiration (excluding the separately recorded BNF carbon cost),
  NPP, LAI, PHU fraction and the leaf/root/mobile/storage carbon stocks.

Stocks are sampled after the daily step; auxiliaries retain the value from their
last daily process stage. The PHU fraction is model development, not an observed
anthesis date. All ordinary harvest/reset days are explicitly marked.

`scripts/attribution/analyze_weather_case.py` is the offline Python analysis.
It checks file hashes, restored-weather masks and the unchanged causal prefix.
It reports weather differences, first resolved process differences, signed peak
and end-of-common-growth stock changes, integrated flux differences over common
growing days, finite yield changes and harvest shifts. Optional static figures
show paired trajectories. Harvest-reset values and missing padding are excluded
from process comparisons; a longer season is not silently compared against zeros.
Photosynthetic capacities have flux-like units but are **not** integrated as
actual carbon fluxes. First-difference thresholds (`rtol=1e-5`, `atol=1e-7`) are
numerical diagnostics, not significance tests or precise causal onset estimates.

This establishes a model-conditional experimental chain: a specified weather
window is changed, process responses are observed, and the finite yield outcome
is measured. It does not identify independent mediated contributions by water,
nitrogen and carbon. Do not sum the one-window finite yield changes or present
them as percentages of the full-reference effect. Interactions and event shifts
prevent that interpretation; moreover full-reference weather remains replaced
beyond the factual harvest, whereas window restorations stop before it.

Before a Paper 1 claim, select observed low-yield events, check that the default
model represents their anomaly, and repeat across explicit reference years and
representative cells. Synthetic validation is not an extreme-event result.

## Reproducible workflow

The isolated runners are in `scripts/attribution`. Operational configuration,
scheduler examples and commands live outside the repository, in the workspace
`workflow/paper1_global_cell/yield_attribution` directory.

The three stages use their own new 600-year HWSD-constrained agricultural
warm-up allocation, default-parameter 1901–2019 history, then selected-cell
attribution. Historical initialisation uses the allocation contract, as in the
existing baseline workflow; it does not claim to restore every final warm-up
state variable. A selected cell is replayed from 1901 and its annual yield must
match the new historical result before its post-sowing state is used for AD.
No old trained profile or production checkpoint is silently substituted.

Each stage records the pinned source, full configuration hash, scientific
configuration hash and upstream allocation/output hash. Case/reference choices
and scheduler resources may change after screening; model inputs, process
settings and source may not. Fresh stage directories are claimed atomically,
and symlinks escaping the dedicated Paper 1 output root are rejected.
