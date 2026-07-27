# Original vs revised Agrocosm — feature parity and a matched simulation

> Compares the standalone LPJmL-derived Agrocosm (base revision `2192dc1f`) against the revised
> Terrarium-based version on a matched ten-year wheat simulation, and records the feature-parity audit.
> The revised model runs with the **full nitrogen cycle wired** (crop demand/uptake kinetics + fertilizer,
> and soil litter-N / immobilization / NH₃ volatilization).

Date: 2026-07-24 (updated 2026-07-27)

## Setup

The same site and forcing are run on both models:

| | Original (standalone) | Revised (Terrarium) |
| --- | --- | --- |
| Crop | `cft1` temperate cereals | `crop_pft("temperate cereals")` |
| Cell | `initial_wheat.jld2` cell 1 | same |
| Forcing | `climate_2000_2009.jld2` (10 yr daily) | same, via `surface_climate_inputs` |
| Initial conditions | `initial_wheat.jld2` (`initialize_simulation`) | `load_crop_initial_conditions` (PHU, sowing, residue, soil C) |
| Soil physics | LPJmL 5-layer water/heat/freeze-thaw | Terrarium soil energy + Richards hydrology |
| Radiation / surface | LPJmL APAR/albedo/PET | Terrarium PALADYN surface energy balance |
| Nitrogen | full demand/uptake + fertilizer (`:yes`) | full demand/uptake kinetics + fertilizer + soil litter-N / immobilization / volatilization |
| Water stress | optimal-λ (cᵢ/cₐ) solver | off (β = 1) — λ solver not yet re-ported |
| Time stepping | discrete daily | continuous, Δt = 600 s |

The original run is `examples/Example_simulation_for_wheat.ipynb`; the revised run is
`examples/wheat_gpp_npp.jl`.

## Matched-simulation results (10 years, annual totals)

| Year | GPP orig | GPP rev | NPP orig | NPP rev | yield rev | peak LAI orig | peak LAI rev |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.16 | 0.23 | 0.10 | 0.17 | 0.012 | 3.21 | 6.98 |
| 2 | 0.54 | 0.28 | 0.32 | 0.21 | 0.010 | 6.97 | 6.98 |
| 3 | 0.36 | 0.26 | 0.20 | 0.19 | 0.012 | 5.05 | 6.97 |
| 4 | 0.50 | 0.20 | 0.31 | 0.15 | 0.012 | 6.97 | 6.97 |
| 5 | 0.87 | 0.36 | 0.53 | 0.27 | 0.013 | 6.97 | 6.98 |
| 6 | 0.44 | 0.21 | 0.27 | 0.15 | 0.011 | 6.97 | 6.97 |
| 7 | 0.36 | 0.28 | 0.18 | 0.20 | 0.012 | 6.14 | 6.98 |
| 8 | 0.37 | 0.22 | 0.21 | 0.16 | 0.011 | 5.31 | 6.98 |
| 9 | 0.49 | 0.33 | 0.26 | 0.24 | 0.012 | 6.75 | 6.97 |
| 10 | 0.41 | 0.20 | 0.25 | 0.14 | 0.012 | 6.97 | 6.97 |
| **mean** | **0.451** | **0.256** | **0.263** | **0.188** | **0.012** | — | — |

(GPP/NPP in kgC/m²/yr; grain yield in kgC/m².)

## Assessment

- **Same order of magnitude and similar dynamics.** Both models produce a temperate-wheat carbon cycle of
  a few tenths of a kgC/m²/yr, peaking in the warmest year (year 5) and dipping in the coldest, with no
  blow-ups over ten years. The harvest index produces a grain crop (~0.012 kgC/m²/yr) each season.
- **Revised ≈ 0.57× the original on GPP, ≈ 0.71× on NPP.** NPP is the fairer cross-model metric:
  Terrarium reports **GPP net of leaf respiration**, whereas the original `gpp` is the gross flux, so a
  definitional offset inflates the GPP gap. On NPP (both *net* primary production) the models agree to
  within ~29 %.
- **LAI differs systematically:** the revised canopy reaches the phenological maximum (≈ 7.0) every year,
  while the original's peak LAI drops to 3–5 in poorer years. The original caps LAI by the running
  carbon/water/nitrogen deficit (the LAI–NPP feedback); the revised uses the heat-unit LAI trajectory
  only. **This is the single biggest remaining structural difference** and explains much of the residual
  gap.

**Component-level agreement holds where the inputs are controlled:** the ported C3/C4 photosynthesis
reproduces Terrarium's `LUEPhotosynthesis` to `rtol = 1e-10` (`test/crop/test_photosynthesis.jl`), and the
soil C–N transforms are unit-tested against their LPJmL forms. The end-to-end gap is in the *coupling*
(LAI feedback, soil water/temperature, radiation), not the crop-physiology kernels.

## The nitrogen cycle (wired 2026-07-27)

The revised model carries the full LPJmL nitrogen cycle — the crop N supply and the soil N transforms it
draws on. All terms are tested for mass conservation.

**Crop side** (`src/crop/nitrogen.jl`, `photosynthesis.jl`):

- **Demand + uptake kinetics.** `CropNitrogenDemand` sizes the crop's demand from the potential
  (light-derived) Vc_max and the organ carbon pools; `CropNitrogenUptakeKinetics` supplies it by
  Michaelis–Menten root uptake from the **soil mineral-N pools over the root zone** (a layer-thickness-
  weighted column integral of `soil_ammonium`/`soil_nitrate`), with a soil-temperature response and a
  root-carbon factor. Uptake fills the demand deficit, capped by the soil supply.
- **Vcmax nitrogen limitation.** `nitrogen_supported_vcmax` caps `Vc_max = min(Vc_max, nitrogen-supported
  capacity)` — the faithful, *light-dependent* LPJmL Rubisco limitation.
- **Fertilizer** (`examples/wheat_gpp_npp.jl`): a continuous mineral-N application (~30 kgN/ha/yr split
  evenly NH₄/NO₃), matching the original's `fertilizer = :yes`.

**Soil side** (`src/crop/soil_biogeochemistry.jl`):

- **Litter nitrogen pool** (`litter_nitrogen` prognostic). Crop litterfall N and harvest residue N enter
  this *organic* pool (not straight to mineral ammonium) and mineralize to NH₄ only as the litter carbon
  decomposes, at the litter's own C:N.
- **Immobilization** (`CropNitrogenMineralization`). When decomposing litter is nitrogen-poor relative to
  the soil C:N, microbial demand *immobilizes* mineral N — a Michaelis–Menten sink on NH₄.
- **NH₃ volatilization** (`CropVolatilization`). An ammonia sink on the **top soil layer** (`k = Nz`),
  driven by air temperature, wind speed, and soil pH.

**Two findings from the wiring:**

1. **The supply and the limitation are coupled.** Wiring the Vcmax *limitation* first, without the N
   *supply*, dropped GPP to 0.18 — the correct light-dependent physics revealed the crop was
   nitrogen-limited under the earlier first-order N closure. Wiring the demand/uptake kinetics **and**
   fertilizer lifted leaf nitrogen above the structural minimum and recovered GPP to 0.256. The lesson is
   that these primitives must be wired together, not piecemeal.
2. **The soil-N transforms don't change this (fertilized) run's carbon.** Once fertilizer lifts the crop
   out of nitrogen limitation, `nitrogen_capacity ≥ potential_vcmax` throughout, so carbon uptake is set
   by light/temperature/phenology and is insensitive to the exact mineral-N budget. The litter pool,
   immobilization, and volatilization correct the **mineral-N budget** — they matter for unfertilized or
   N-stressed runs, the N balance, and downstream leaching — not this run's GPP/NPP.

## Sources of the residual difference

In rough order of contribution:

1. **Missing λ water-coupling solver (crop physiology — closeable).** The original couples photosynthesis
   to soil water through the optimal-λ (cᵢ/cₐ) bisection (`solve_lambda_*` + `lpj_bisect`); the revised
   uses a crude `λ = λ_min + (λ_opt − λ_min)·β`, and **water stress is off by default** (β = 1,
   `plant_available_water = nothing`, default `SoilHydrology` = `NoFlow`). This severs the
   water→GPP/transpiration feedback central to LPJmL — the single largest remaining value-changing gap.
2. **LAI–NPP carbon-deficit feedback (deferred).** The original caps LAI by the running carbon/water/N
   deficit; the revised follows the heat-unit LAI trajectory only, so its canopy sits at the phenological
   maximum every year. This drives much of the LAI (and hence light-absorption) difference.
3. **Vernalization + climate buffer/spinup (missing).** Winter wheat needs vernalization for correct
   phenology timing; both were removed and not re-added.
4. **Definitional GPP offset (not a defect).** Terrarium's `gpp` is net of leaf respiration; the original's
   is gross. On NPP (both *net*) the models already agree to ~29 %.
5. **Structural, intended.** A different soil water/temperature model (Richards + energy vs the LPJmL
   5-layer scheme) and surface-radiation scheme (PALADYN surface energy balance vs LPJmL APAR/albedo/PET)
   necessarily change water stress, light absorption, and timing. Exact numerical agreement is not expected
   once the soil/surface is Terrarium's.

## Feature-parity audit (base revision `2192dc1f` vs current tree)

Status legend: **PRESENT** (wired into the running model), **REPLACED-BY-TERRARIUM**, **MISSING/DEFERRED**.

**Crop physiology:** C3/C4 photosynthesis PRESENT · APAR PRESENT · carbon allocation PRESENT ·
autotrophic respiration PRESENT · phenology + LAI PRESENT · crop carbon PRESENT · crop nitrogen PRESENT
(full demand/uptake kinetics) · Vcmax N-limitation PRESENT (faithful, light-dependent) · harvest index
PRESENT · temperature stress PRESENT · root distribution PRESENT · **λ solver MISSING** ·
**vernalization MISSING** · **LAI↔NPP carbon-deficit feedback DEFERRED**.

**Management:** sowing PRESENT · harvest PRESENT (harvest-index yield) · fertilizer PRESENT (continuous
flux) · **tillage MISSING/DEFERRED**.

**Soil biogeochemistry:** decomposition PRESENT · nitrification PRESENT · denitrification PRESENT ·
mineralization PRESENT · litter routing PRESENT · litter-nitrogen pool PRESENT · immobilization PRESENT ·
NH₃ volatilization PRESENT (top-layer sink) · **NO₃ leaching MISSING** · **surface-litter
hydrology/thermal MISSING**.

**Soil/surface (delegated to Terrarium):** soil temperature, evaporation, infiltration/percolation,
transpiration, interception, freeze-thaw, radiation/albedo/PET all REPLACED-BY-TERRARIUM (note: default
crop model uses `NoFlow` hydrology → static water) · **snow MISSING** (Terrarium has only an abstract
stub).

**Climate/forcing:** climate input REPLACED-BY-TERRARIUM (`surface_climate_inputs` → `FieldTimeSeries`) ·
CO₂ PRESENT (constant; time-series DEFERRED) · **climate buffer/spinup MISSING**.

**Diagnostics:** **runtime Water/Nitrogen/Carbon/Thermal balance ledgers MISSING** (conservation is
unit-tested only).

**Numerics:** CPU PRESENT · GPU PRESENT (framework; full-model run DEFERRED) · Float32/Float64 PRESENT ·
checkpoints REPLACED-BY-TERRARIUM · differentiability PARTIAL (soil biogeochem via Reactant; full crop
`LandModel` blocked by the root-fraction Reactant gap).

### Remaining gaps, ordered by impact on a wheat run

1. λ water-coupling solver — MISSING (water→GPP feedback severed; water stress off by default).
2. LAI↔NPP carbon-deficit feedback — DEFERRED (revised canopy sits at the phenological maximum).
3. Vernalization + climate buffer/spinup — MISSING (winter-wheat phenology timing).
4. NO₃ leaching — MISSING.
5. Surface-litter hydrology/thermal — MISSING.
6. Tillage — MISSING (upstream Terrarium hook).
7. Snow model — MISSING (Terrarium stub only).
8. Runtime balance-ledger diagnostics — MISSING.

The λ solver (#1), the LAI–NPP feedback (#2), and vernalization (#3) are the remaining **value-changing**
gaps; all three require re-porting removed code rather than wiring a tested primitive. #4–#8 are narrower
fidelity/robustness gaps. These are tracked in `2026-07-24_NOTES_future_work.md`.

## Reproduce

```
julia --project=. examples/wheat_gpp_npp.jl        # revised, 10-year wheat GPP/NPP/yield
```

The original numbers are from `examples/Example_simulation_for_wheat.ipynb` on base revision `2192dc1f`.
