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
- **Open:** a concrete source-order or formulation difference still needs a
  decision before claiming LPJmL-style behavior.

## Daily process order

The two drivers are `src/simulations/daily_crop_C3.jl` and
`src/simulations/daily_crop_C4.jl`.  Their common scientific pathway is:

```text
climate history → cultivation/tillage/bioturbation → albedo/PET → snow
     → soil hydraulic/litter/thermal preparation
     → pre-crop soil C–N decomposition → phenology/harvest/residue routing
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
| Stand termination | harvested GPU sentinel reconstructed by `cultivate!` | `killstand()`, `delpft()` | Adapted | LPJmL deletes a crop PFT; Agrocosm zeros/reconstructs seasonal state in a fixed allocation. The resulting inactive-crop behavior is tested, while the representation is intentionally GPU-oriented. |
| Snow, albedo, and PET | `albedo!`, `petpar!`, `snow!` | `albedo_stand()`, `albedo_crop()`, `petpar()`, `snow.c` | Aligned / Adapted | Agrocosm preserves LPJmL's albedo/PET-before-snow order and now reconstructs the full green-canopy, surface-litter, bare-soil, and start-of-day snow mixture. Inactive crop cells use the bare stand soil/snow mixture. The fixed-array kernel and direct reconstruction of litter cover from current carbon are GPU-oriented adaptations that avoid stale cached cover. Current snow still feeds same-day soil thermal resistance and canopy radiation later in the step. |
| Tillage and bioturbation | `litter_tillage!`, `litter_bioturbation!` | `cultivate.c`, `update_daily_cell.c` | Aligned | Cultivation-tillage is event based and the daily surface-to-subsurface transfer precedes soil decomposition. |
| Soil physical preparation | `pedotransfer!`, `update_surface_litter_properties!`, `soil_temperature!` | soil thermal update, `pedotransfer()`, `updatelitterproperties()` | Aligned / Adapted | The hydraulic/litter pair now follows LPJmL's pedotransfer-before-litter order. Agrocosm intentionally performs both before its enthalpy solver because that solver consumes current pore volume and current litter depth/water as thermal properties; LPJmL's thermal solver is ordered earlier. All remain before decomposition. |
| C–N decomposition | `soil_cn_decomposition!` | `daily_littersom.c`, `littersom_nomethane.c` | Aligned / Adapted | Same pre-crop role: decomposition, respiration, mineralization/immobilization, and nitrification make mineral N available to the crop. Agrocosm uses one annual-crop litter class and shared post-spin-up `c_shift` profiles instead of LPJmL's PFT litter list. |
| Phenology and normal harvest | `phenology_crop!`, `harvest_crop!` | `phenology_crop.c`, `harvest_crop.c` | Aligned | Both run before infiltration and daily crop assimilation. New residues are routed after the day's decomposition, so they begin decomposing the next day. |
| Interception and infiltration | `interception!`, `soil_infiltration!` | `interception.c`, `infil_perc.c` | Aligned / Adapted | Same placement before water-stressed assimilation. Agrocosm retains an explicit five-layer enthalpy ledger and GPU-safe thermal update schedule. |
| C3/C4 assimilation and water limitation | `photosynthesis_C3!` / `photosynthesis_C4!`, `transpiration!`, `solve_lambda_*` | `photosynthesis.c`, crop water-stress path | Adapted | The LPJmL-informed potential-capacity → water-limited-λ → final-assimilation sequence is retained. The GPU-compatible λ solver and numerical guards are an implementation adaptation. |
| Carbon allocation and respiration | `crop_carbon!` | `npp_crop.c`, `allocation_daily_crop.c` | Aligned / Adapted | Agrocosm follows LPJmL's default `crop_resp_fix=true` configuration by using fixed organ N:C ratios. It intentionally deducts root respiration from NPP to close crop carbon, correcting the omission in LPJmL's `npp_crop.c`; organ allocation otherwise follows the LPJmL crop path with a fixed-array representation. |
| Crop N demand, uptake, and allocation | `crop_nitrogen!`, `allocate_crop_nitrogen!` | `ndemand_crop.c`, `nuptake_crop.c`, `vmaxlimit_crop.c` | Adapted | Soil-N supply, uptake, allocation, and management inputs are active. The optional N-to-`vcmax` feedback remains disabled by default until the full LPJmL feedback sequence is completed. |
| Soil evaporation and plant water removal | `evaporation!`, `soil_evapotranspiration!` | `waterbalance.c` | Aligned / Adapted | Both apply after crop demand/growth calculations. Agrocosm keeps explicit daily water-flux arrays for conservation and GPU execution. |
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

## Current boundaries

The audit covers the prescribed, single-crop, non-methane path.  It does not
claim equivalence for multiple concurrent stands, stand-fraction competition,
dynamic sowing, rotation, irrigation infrastructure, wetlands/rice methane,
or global production configuration.  These boundaries are part of the
[roadmap](roadmap.md), rather than hidden approximations.

## Closed daily-order decisions

The C3 and C4 drivers now share the following tested ordering constraints:

1. `update_climbuf!` precedes `cultivate!`.
2. `albedo!` and `petpar!` precede `snow!`.
3. `pedotransfer!` precedes `update_surface_litter_properties!`, which precedes
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
