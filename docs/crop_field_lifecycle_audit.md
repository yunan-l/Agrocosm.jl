# Crop field lifecycle audit

This audit fixes and records the lifecycle contract of the reorganized `Crop`
storage. It covers the normal growing day, fallow day, sowing day, harvest day,
the day after harvest, output, and checkpoint boundaries.

## Daily ownership

Agrocosm does not launch a global daily reset kernel. Every daily field has one
owner process, and that process writes both its active and inactive branches.
This avoids an additional full-grid GPU traversal and keeps lifecycle semantics
next to the scientific operation that produces the value.

| Category | Boundary rule |
|---|---|
| `state` | Preserved. Contains every quantity read by a later daily transition that cannot be reconstructed solely from current forcing, parameters, and static configuration. |
| `fluxes` | Its owner overwrites it every day and writes zero when inactive. |
| `events` | Cultivation and harvest each write either zero or one every day. |
| `auxiliary` | Contains same-day derived values or explicitly named cross-day process memory. |
| `workspace` | Its owner overwrites it before every read. |

The checkpoint boundary is `crop_restart_payload(...).state` plus explicit
static process configuration. Implementation workspace is never included.

## Corrected classifications

The audit found the following earlier classification errors:

1. `lai`, `laimax_adjusted`, and `lai_npp_deficit` are
   `state.canopy`: next-day interception, radiation, or LAI growth reads them
   before they can be recomputed. `phenology_fraction = lai / laimax` is a
   scalar local workspace inside the albedo kernel. The LPJmL-compatible nonnegative
   `actual_lai = max(0, lai - lai_npp_deficit)` remains a daily auxiliary.
2. N and water sufficiency are now `state.nitrogen.sufficiency` and
   `state.water.sufficiency`, because the next day's LAI update consumes them.
   Daily demands and deficit factors remain `auxiliary.stress`.
3. `fphu = husum / phu` is a daily diagnostic. Fertilizer and phenology derive
   it directly from the prognostic heat sum rather than trusting a prior cache.
   `phu` and `winter_type` are static process configuration.
4. `sowing_date` is prescribed configuration. Harvest date, annual yield, and
   annual-harvest presence are output bookkeeping, not crop state.
5. `output.annual` holds yield and harvest date until the day-365 annual output
   row is written; daily `fluxes.carbon.yield` remains zero off harvest day.
6. `fluxes.water.transpiration` was never written. It was removed; the complete
   water flux is `transpiration_layer`, whose layer sum supplies diagnostics.

`fluxes.carbon.harvest_export` was added so carbon conservation does not depend
on reading plant stocks after harvest-day owner processes clear them.

## Field ownership

| Group | Fields | Owner/overwrite stage |
|---|---|---|
| carbon flux | `yield`, `harvest_export` | harvest |
| carbon flux | `gross_assimilation`, `net_assimilation`, `water_limited_assimilation`, `leaf_respiration` | photosynthesis |
| carbon flux | `respiration` | respiration |
| carbon flux | `npp` | carbon allocation |
| nitrogen flux | `seed_input` | cultivation |
| nitrogen flux | prescribed manure/fertilizer | fertilizer |
| nitrogen flux | `uptake`, `auto_fertilizer` | nitrogen uptake |
| nitrogen flux | `harvest_export` | harvest |
| water flux | `interception` | interception |
| water flux | `transpiration_layer` | transpiration |
| event | `sowing` | cultivation, exactly on `sowing_date` |
| event | `harvest` | harvest transition `!harvesting_previous && harvesting` |
| phenology state | `vdsum`, `husum`, phase flags and growing-day count | cultivation and phenology |
| canopy state | `lai`, `laimax_adjusted`, `lai_npp_deficit` | cultivation, phenology, LAI allocation |
| phenology diagnostic | `fphu = husum / phu` | phenology, fertilizer, carbon allocation, output |
| phenology/calendar configuration | `phu`, `winter_type`, `sowing_date` | initialization and cultivation |
| output bookkeeping | annual yield and harvest date; annual-harvest flag is derived at reporting | harvest and day-365 output |
| canopy auxiliary | current-day `flaimax`, actual LAI, albedo, radiation and conductance | phenology, albedo, radiation, interception, transpiration |
| photosynthesis auxiliary | capacity, lambda and temperature limitation | temperature stress, photosynthesis and lambda solve |
| stress/root auxiliary | same-day demands and deficit factors; root distribution and top-three-layer root-zone available water | nitrogen and water diagnostics |
| workspace | currently empty | reserved for future preallocated kernel scratch arrays |

The owner-overwrite invariant is tested by poisoning every daily flux,
recomputable auxiliary field, event and workspace between two simulation
blocks. Explicit cross-day process memory is excluded from poisoning. The
second day must reproduce an unpoisoned simulation exactly.

## Harvest equivalence to LPJmL stand deletion

The harvest sequence is now explicit:

1. detect the transition and set `events.harvest`;
2. record yield, total carbon export, nitrogen export and litter inputs;
3. add grain yield to the stand-independent annual calendar accumulator;
4. preserve the established harvest-day diagnostic/output values;
5. let the existing carbon, nitrogen, water, canopy and phenology owner kernels
   write their inactive values; no additional retirement kernel is launched.

Live carbon and nitrogen stocks and daily fluxes are cleared on the harvest day
by their normal inactive branches. Canopy and phenology memory is cleared on
the following daily transition. Multiplicative N and water process-memory values become
one (neutral). Crop configuration remains: `phu`, `winter_type`, and the
sowing calendar. Harvest date and annual yield reside separately in output
bookkeeping and never feed the crop transition.
The inactive crop placeholder is therefore scientifically equivalent to an
absent LPJmL stand without paying for an extra full-grid kernel every day.

## Output and restart

Output remains an explicit selection of scientific state, flux and auxiliary
fields. No workspace field is present in `Output` or `CropOutput`.

`crop_restart_payload(crop)` exposes all checkpointed crop state plus static
phenology/calendar configuration. Daily fluxes, events, recomputable
diagnostics and workspace are excluded. Root distribution is reconstructed
from PFT parameters at initialization. A full `CropSimulation` checkpoint must
also retain `output.annual` when it may resume before the year-end annual row
is emitted. A restart must be written only at a completed daily boundary.

## Regression coverage

A 730-day Float32 reference was generated directly from pre-refactor commit
`f6059d5`. The refactored model was compared against all 146 arrays in `Output`
and final `Soil` storage. The maximum absolute difference was `0.0` (bitwise
identical on the reference CPU).

An additional comparison stops on active growing day 120 and maps fields by
scientific meaning across the old and new layouts. All mapped crop, soil and
output arrays were bitwise identical before the subsequent scientific
corrections described below.

The lifecycle audit also removed the separate daily harvested-stand retirement
kernel. Existing inactive branches own state cleanup, and the harvest flags are
cleared inside the already-running phenology kernel. Root-zone water output was
fused into the already-running harvest kernel, removing another daily kernel
launch without changing its three-layer calculation or timing.

After establishing that bitwise structural baseline, a separate scientific
correction aligned LAI-deficit handling with LPJmL. Potential phenological
`state.canopy.lai` is no longer reduced in place each day. Consumers and output
instead use nonnegative `auxiliary.canopy.actual_lai =
max(0, lai - lai_npp_deficit)`. This deliberately changes trajectories affected
by an NPP deficit and prevents the former multi-day repeated subtraction that
produced negative LAI; it is not a consequence of the storage refactor.

Further scientific corrections align LPJmL temperature stress, respiration
temperature sources, nonnegative growth respiration and snow-cover suppression
of canopy absorbed PAR. Root respiration is also subtracted from NPP, correcting
the current LPJmL crop carbon bookkeeping after consultation with its developers.
These changes are independently tested
against LPJmL formulas and can intentionally change trajectories; they are not
storage-refactor differences. Negative-biomass crop termination remains
deferred by design.

- `test_field_lifecycle.jl`: CPU Float32/Float64 poisoning test proving that
  owner processes overwrite daily fields without a global reset;
- `test_checkpoint_minimal_state.jl` and its GPU counterpart: verifies that
  only process-changing fields enter `crop.state`, `fphu` is reconstructed
  from `husum / phu`, albedo derives its local canopy fraction from `lai`, and
  annual harvest records remain output bookkeeping;
- `test_crop_lifecycle.jl`: two complete crop seasons, event uniqueness,
  fallow outputs and post-harvest state;
- `test_harvest_balance.jl`: conservative C/N export and residue routing;
- `test_crop_lifecycle_gpu.jl`: CPU/GPU lifecycle equivalence.
- `test_actual_lai.jl` and `test_actual_lai_gpu.jl`: potential-LAI persistence
  and nonnegative actual LAI in carbon allocation, radiation and interception.
- `test_temperature_stress_lpjml.jl`: independent C3/C4 LPJmL temperature
  response formula.
- `test_respiration_lpjml.jl`: air/soil temperature sources, nonnegative growth
  respiration, and corrected NPP subtraction of complete respiration including
  root respiration.
- `test_canopy_snow_cover.jl` and its GPU counterpart: snow-covered crop canopy
  has zero absorbed PAR for both C3 and maize paths.
