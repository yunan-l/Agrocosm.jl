# LPJmL process and daily-order audit

This is a source-level audit of Agrocosm's current single-crop pathway against
LPJmL commit `572e2b906ac2c55b2ee6661a93e4633b126254e4`.  It is deliberately
not an assertion that every LPJmL internal variable or every floating-point
operation is reproduced.  The purpose is to make the scientific lineage,
daily data flow, and intentional simplifications explicit.

## Status labels

- **Aligned:** the process and its position in the daily sequence agree at the
  level relevant to the current single-crop model.
- **Adapted:** the same process is present, but its representation is
  intentionally different for the current CPU/GPU architecture.
- **Corrected:** a source comparison found a behavior that did not match the
  intended LPJmL rule; the implementation and a focused regression test were
  corrected together.
- **Open:** a concrete source-order or formulation difference still needs a
  decision before claiming LPJmL-style behavior.

## Daily process order

The shared C3/C4 driver is `src/simulations/daily_crop.jl`. Compile-time
pathway dispatch selects the C3 or C4 radiation, photosynthesis, and water-
limited lambda operations without changing their common scientific sequence:

```text
climate history → cultivation/tillage/bioturbation → albedo/PET → snow
     → soil hydraulic/litter/thermal preparation
     → pre-crop soil C–N decomposition → pre-phenology raw gp snapshot (N path)
     → phenology/harvest/residue routing
     → interception/infiltration → photosynthesis and water limitation
     → crop carbon and nitrogen → evaporation/water removal
     → denitrification and NH3 volatilization
```

LPJmL splits the equivalent work between `update_daily_cell.c`,
`daily_littersom.c`, and `daily_agriculture.c`: it updates the climate buffer,
sows/kills stands, applies bioturbation and PET, processes snow, performs soil
thermal and litter/SOM work, then calls `daily_stand()` (and therefore
`daily_agriculture()`), followed by denitrification and volatilization.

| Stage | Agrocosm entry point | LPJmL source basis | Status | Audit conclusion |
| --- | --- | --- | --- | --- |
| Climate history and sowing | `update_climbuf!`, `cultivate!` | `daily_climbuf()`, `sowing()` | Aligned | Agrocosm now updates the climate buffer before cultivation. Current prescribed sowing does not consume that update, but the dependency direction is locked by the daily-order contract test for future dynamic sowing. |
| Stand termination | harvested GPU sentinel reconstructed by `cultivate!` | `killstand()`, `delcft()` | Adapted | LPJmL deletes a crop CFT; Agrocosm zeros/reconstructs seasonal state in a fixed allocation. The resulting inactive-crop behavior is tested, while the representation is intentionally GPU-oriented. |
| Snow, albedo, and PET | `albedo!`, `petpar!`, `snow!` | `albedo_stand()`, `albedo_crop()`, `petpar()`, `snow.c` | Aligned / Adapted / Corrected | Agrocosm preserves LPJmL's albedo/PET-before-snow order and reconstructs the full green-canopy, surface-litter, bare-soil, and start-of-day snow mixture. Inactive crop cells use the bare stand soil/snow mixture. PET now retains LPJmL's lower-zero clamp without the former non-source 15 mm d⁻¹ upper cap. The fixed-array kernel and direct reconstruction of litter cover from current carbon are GPU-oriented adaptations that avoid stale cached cover. Current snow still feeds same-day soil thermal resistance and canopy radiation later in the step. |
| Tillage and bioturbation | `litter_tillage!`, `tillage_hydraulics!`, `litter_bioturbation!` | `cultivate.c`, `tillage.c`, `update_daily_cell.c` | Aligned / Adapted | Cultivation moves litter and reduces the top-layer bulk-density factor; accepted infiltration subsequently settles that factor toward one. Agrocosm's current single-crop path assumes tillage is enabled and omits LPJmL's water-table gate because it has no prognostic water table. |
| Soil physical preparation | `pedotransfer!`, `update_surface_litter_properties!`, `soil_temperature!` | soil thermal update, `pedotransfer()`, `updatelitterproperties()` | Aligned / Adapted | `pedotransfer!` now applies LPJmL's tillage correction to top-layer saturation, field capacity, retention exponent, holding capacity, and saturated conductivity while conserving absolute water and ice stocks. The hydraulic/litter pair follows LPJmL's pedotransfer-before-litter order. Agrocosm intentionally performs both before its enthalpy solver because that solver consumes current pore volume and current litter depth/water as thermal properties; LPJmL's thermal solver is ordered earlier. |
| C–N decomposition | `soil_cn_decomposition!` | `daily_littersom.c`, `littersom_nomethane.c` | Aligned / Adapted | Same pre-crop role: decomposition, respiration, mineralization/immobilization, and nitrification make mineral N available to the crop. Agrocosm uses one annual-crop litter class and shared post-spin-up `c_shift` profiles instead of LPJmL's CFT litter list. |
| Phenology and normal harvest | `phenology_crop!`, `harvest_crop!` | `phenology_crop.c`, `harvest_crop.c` | Aligned / Corrected | Both run before infiltration and daily crop assimilation. Actual daily LAI increment now uses LPJmL's `min(wscal, vscal)` exactly; the former additional `/ 1.5` water scaling was removed. New residues are routed after the day's decomposition, so they begin decomposing the next day. |
| Interception and infiltration | `interception!`, `soil_infiltration!` | `interception.c`, `infil_perc.c` | Aligned / Adapted | Same placement before water-stressed assimilation. The bounded 4 mm slug loop now uses LPJmL's `MAXITER = 1000`; if it reaches the cap, remaining infiltration is routed to surface runoff rather than left unaccounted. Agrocosm retains an explicit five-layer enthalpy ledger and GPU-safe thermal update schedule. |
| C3/C4 assimilation and water limitation | `photosynthesis_C3!` / `photosynthesis_C4!`, `transpiration!`, `solve_lambda_*` | `photosynthesis.c`, `gp_sum.c`, crop water-stress path | Adapted / Corrected | The N-aligned path snapshots raw `gp_sum` from the pre-phenology canopy and uses it for the initial water pass, while the current-canopy pass still supplies the λ solve. No additional `gp × fpar/phen` scaling is applied: it is absent from 5.10 and duplicates canopy scaling already present in APAR/conductance. The GPU-compatible λ solver and numerical guards remain implementation adaptations. |
| Carbon allocation and respiration | `crop_carbon!` | `npp_crop.c`, `allocation_daily_crop.c` | Aligned / Adapted | Agrocosm defaults to the released LPJmL 6.1.9 configuration's `crop_resp_fix=true`, using fixed organ N:C for crop respiration; the dynamic organ-N:C mode remains an explicit opt-out. The 6.1.9 dynamic storage/pool cap assignments are internally inconsistent and the cap block is absent from 5.7.9, so Agrocosm preserves its existing opt-out behavior rather than copying that extreme-state branch. NPP is defined after leaf, root, storage, pool, and growth respiration; organ allocation uses the LPJmL crop pathway as a historical process basis with a fixed-array representation. |
| Crop N demand, uptake, and allocation | `crop_nitrogen!`, `allocate_crop_nitrogen!` | `ndemand_crop.c`, `nuptake_crop.c`, `vmaxlimit_crop.c` | Adapted | Soil-N supply, uptake, allocation, and management inputs are active. The optional N-to-`vcmax` feedback remains disabled by default until the full LPJmL feedback sequence is completed. |
| Soil evaporation and plant water removal | `evaporation!`, `soil_evapotranspiration!` | `waterbalance.c` | Aligned / Adapted / Corrected | Both apply after crop demand/growth calculations. Managed land now uses 6.1.9's `0.1` minimum exposed fraction in both evaporation-energy and litter-cover factors. As in 5.10/6.1, liquid plus ice drives the sigmoid moisture response, while soil evaporation is capped and distributed using only liquid water remaining after transpiration. The low-level explicit legacy switch is retained only for regression tests. Agrocosm keeps explicit daily water-flux arrays for conservation and GPU execution. |
| Late N losses | `post_crop_nitrogen_losses!` | `denitrification.c`, `volatilization.c` | Aligned | Both occur after the daily stand/crop update, using the updated mineral-N pools and moisture state. |

## Process coverage by domain

### Crop

| Process family | Agrocosm sources | LPJmL scientific basis | Status |
| --- | --- | --- | --- |
| C3/C4 photosynthesis, temperature stress, APAR, λ water limitation | `processes/crop/photosynthesis.jl`, `transpiration.jl`, `lambda_solver.jl`, `radiation.jl` | `photosynthesis.c`, crop water-stress pathway | Adapted |
| Phenology, canopy, sowing and harvest | `phenology.jl`, `cultivate.jl`, `harvesting.jl`, `albedo.jl` | `phenology_crop.c`, `cultivate.c`, `harvest_crop.c`, `albedo_crop.c` | Aligned / Adapted |
| Carbon growth, respiration, and allocation | `crop_carbon.jl`, `respiration.jl`, `carbon_allocation.jl` | `npp_crop.c`, `allocation_daily_crop.c` | Aligned / Adapted |
| Crop nitrogen and fertilization | `nitrogen_allocation.jl`, `fertilizer.jl` | `ndemand_crop.c`, `nuptake_crop.c`, `vmaxlimit_crop.c` | Adapted |

### Soil

| Process family | Agrocosm sources | LPJmL scientific basis | Status |
| --- | --- | --- | --- |
| Snow, multilayer water, freezing/thawing, and energy transport | `processes/climate/snow.jl`, `processes/soil/soil_water.jl`, `soil_temperature.jl` | `snow.c`, `infil_perc.c`, soil thermal routines | Adapted |
| Litter routing and soil C–N decomposition | `litter_routing.jl`, `soil_carbon.jl`, `soil_nitrogen.jl` | `daily_littersom.c`, `littersom_nomethane.c` | Aligned / Adapted |
| N transformations and losses | `nitrogen_transform.jl` | `denitrification.c`, `volatilization.c` | Aligned |

## 2026-07 source comparison pass

This pass compared the active daily crop route against the named LPJmL source
files, rather than only checking numerical outputs.  It covers the processes
that execute for the present prescribed annual-crop configuration.

| Process group | LPJmL source inspected | Conclusion |
| --- | --- | --- |
| Climate history, sowing, phenology, harvest, residue timing | `daily_climbuf()`, `cultivate.c`, `phenology_crop.c`, `harvest_crop.c` | Daily order and prescribed-crop behavior align. Dynamic sowing and crop rotations remain outside the present mode. |
| Albedo, PET, snow, interception | `albedo_crop.c`, `petpar.c`, `snow.c`, `interception.c` | Snow phase partition, fixed sublimation, melt, and geometry follow the LPJmL pathway. The PET upper cap was a real discrepancy and is corrected below. |
| Photosynthesis, water stress, lambda solve | `photosynthesis.c`, `gp_sum.c`, `water_stressed.c` | Potential assimilation, the pre-phenology raw-conductance snapshot, water-limited lambda closure, and the 30-step primary bisection are consistent with the version-screened crop path. The extra 6.1.9 `gp × fpar/phen` operation is deliberately rejected because it is absent in 5.10 and duplicates existing canopy scaling. Water-table/inundation terms are not represented because Agrocosm has no prognostic water table. |
| Leaf development, respiration, carbon allocation and crop failure | `phenology_crop.c`, `npp_crop.c`, `allocation_daily_crop.c` | LAI water scaling is corrected below. Allocation, senescence pool draw-down, and negative-biomass termination follow the LPJmL control flow. Agrocosm defines NPP after all crop respiration components. |
| N demand, uptake, allocation, mineralization, nitrification, denitrification, NH₃ | `ndemand_crop.c`, `nuptake_crop.c`, `daily_littersom.c`, `littersom_nomethane.c`, `denitrification.c`, `volatilization.c` | Current fixed-N:C crop-respiration mode and daily ordering are consistent. N-to-`vcmax` feedback stays deliberately disabled by default. Ice treatment in denitrification is an explicit frozen-soil simplification. |
| Tillage, pedotransfer, infiltration/percolation, evaporation | `cultivate.c`, `pedotransfer.c`, `infil_perc.c`, `waterbalance.c` | Hydraulic parameter update and daily placement are retained. The infiltration cap was a real discrepancy and is corrected below. Full irrigation remains the configured LPJmL-style field-capacity reset, not a water-supply model. |
| Soil thermal / freeze-thaw transport | LPJmL soil thermal routines | Not claimed equivalent: frozen-soil infiltration and heat transport are deliberately deferred from the current scope. The existing enthalpy solver is conservation-tested independently. |

### Corrections found by the source pass

| Process | Former behavior | LPJmL-compatible behavior now used | Regression evidence |
| --- | --- | --- | --- |
| Infiltration/percolation | Stopped after 500 4-mm slugs and could leave the remaining daily input unaccounted. | Uses `MAXITER = 1000`, then routes any remainder to surface runoff and zeros the mutable infiltration input, as in `infil_perc.c`. | `test/processes/soil/test_soil_water.jl` checks the cap and daily water accounting. |
| Actual LAI development | Multiplied potential LAI growth by `min(wscal / 1.5, vscal)`. | Uses `min(wscal, vscal)`, matching `phenology_crop.c`; both factors are already 0–1. | `test/processes/crop/test_actual_lai.jl` checks a water-limited canopy increment. |
| Equilibrium evaporation | Clamped PET to `[0, 15]` mm d⁻¹. | Keeps LPJmL's lower-zero clamp only; high-energy values are not artificially capped. | `test/processes/crop/test_radiation_lpjml.jl` checks a value above 15 mm d⁻¹. |
| Crop conductance timing | Recomputed the initial water-pass conductance from the canopy after the day's phenology update. | The N-aligned path preserves LPJmL 5.10/6.1's pre-phenology raw `gp_sum`, without 6.1.9's additional `fpar/phen` rescaling. | `test/processes/crop/test_lambda_water_coupling.jl` changes the post-phenology canopy inputs and verifies that the stored raw conductance drives demand and water sufficiency. |
| Managed-land soil evaporation | Used a `0.05` minimum exposed fraction and could diagnose more soil evaporation than the liquid water remaining after transpiration. | Uses 6.1.9's `0.1` managed-land floor in both exposed-fraction terms; liquid plus ice supplies `w_evap`, while post-transpiration liquid water supplies the cap and ratio denominator (`tmpwater`). | `test/processes/soil/test_surface_litter_water.jl` covers both floors, frozen-water response, the dry-soil cap, and the explicit legacy branch. |

The scientific temperature coefficient remains `k_temp = 0.0693`, as declared
in the LPJmL parameter file and read by 5.10. LPJmL 6.1.9's scanner omits this
field and therefore leaves its global static value at zero; Agrocosm treats
that executable behavior as a reader bug rather than a parameter target.

Nitrate transport intentionally retains the 5.10 timing: NO₃ moves within
each bounded `≤ 4 mm` infiltration substep. Although 6.1.9 defers transport to
one post-hydrology pass, it resets the local lateral-runoff accumulator inside
each substep before using that value at day end, so that path is not adopted as
a scientific reference. Agrocosm therefore keeps the internally consistent
5.10 operator and its paired `NPERCO = 0.4`, rather than mixing that time scale
with 6.1.9's `0.3`. This decision is independent of the `NO_METHANE` soil C–N
order.

Nitrification moisture coefficients are soil properties, not global model
parameters. Agrocosm maps the LPJ soil classes SaCl, SaClLo, SaLo, LoSa, and Sa
(`soilcode` 3, 6, 9, 11, and 12) to `(0.55, 1.7, -0.007, 3.22)` and all other
classes to `(0.45, 1.27, 0.0012, 2.84)`, matching `soil.cjson` and `f_wfps()`.
The mapping uses the discrete LPJ `soilcode` already carried by the production
input contract. It does not infer an LPJ class from continuous HWSD texture;
inputs without a soil code retain the fine/default tuple for compatibility.
The same initialization maps LPJ's anion-excluded porosity to `0.4` for those
five sandy classes plus clay-light/rock, and to `0.3` otherwise; nitrate
transport reads that per-cell property while retaining the selected 5.10
substep timing and `NPERCO = 0.4`.

The production soil NetCDF maps source code 13 to rock/ice, whereas Agrocosm's
canonical 14-class lookup uses code 13 for clay-light and 14 for rock. The
AgrocosmData reader therefore resolves the NetCDF `map` names into canonical
codes instead of assuming that source numbers are lookup positions. The two
classes share the N coefficients above, so this correction does not alter the
N-parameter mapping. It does alter their hydraulic and thermal properties.
The ten-grid fixture and cell 60866 contain no source-code-13 cells, but the
global wheat inputs do; derived initial-data and checkpoint artifacts that
were built with the old positional interpretation must be regenerated.

The infiltration loop retains Agrocosm's `1e-5` residual threshold rather than
LPJmL's smaller floating-point epsilon. This is a numerical tolerance choice;
after the cap fallback it cannot discard a water input and is not treated as a
separate process discrepancy.

## Current boundaries

The audit covers the prescribed, single-crop, non-methane path.  It does not
claim equivalence for multiple concurrent stands, stand-fraction competition,
dynamic sowing, rotation, irrigation infrastructure, wetlands/rice methane,
or global production configuration.  These boundaries are part of the
[roadmap](roadmap.md), rather than hidden approximations.

## Closed daily-order decisions

The C3 and C4 drivers now share the following tested ordering constraints:

1. `update_climbuf!` precedes `cultivate!`.
2. `litter_tillage!` and `tillage_hydraulics!` precede `pedotransfer!`.
3. `albedo!` and `petpar!` precede `snow!`.
4. `pedotransfer!` precedes `update_surface_litter_properties!`, which precedes
   `soil_temperature!`.

The lightweight `test/simulations/test_daily_process_order.jl` contract test
guards orchestration order. Existing process tests continue to cover snow,
surface-litter water conservation, hydraulic repartitioning, and snow/litter
thermal resistance numerically.

## Next audit actions

1. Once a public checkpoint restore API exists, test uninterrupted versus
   save/restore continuation on CPU and CUDA, rather than only payload shape.
2. Re-run this audit whenever a new selectable process model is introduced in
   Phase 2.
