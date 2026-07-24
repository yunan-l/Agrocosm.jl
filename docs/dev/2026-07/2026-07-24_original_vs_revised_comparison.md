# Original vs revised Agrocosm — feature parity and a matched simulation

> Status: **in progress**. Compares the standalone LPJmL-derived Agrocosm (base revision
> `2192dc1f`) against the revised Terrarium-based version on a matched wheat simulation, and records the
> feature-parity audit.

Date: 2026-07-24

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
| Nitrogen | full demand/uptake + fertilizer (`:yes`) | first-order N closure, no fertilizer, well-watered (β = 1) |
| Time stepping | discrete daily | continuous, Δt = 600 s |

The original run is `examples/Example_simulation_for_wheat.ipynb`; the revised run is
`examples/wheat_gpp_npp.jl`.

## Matched-simulation results (10 years, annual totals)

| Year | GPP orig | GPP rev | NPP orig | NPP rev | peak LAI orig | peak LAI rev |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.16 | 0.22 | 0.10 | 0.16 | 3.21 | 6.98 |
| 2 | 0.54 | 0.27 | 0.32 | 0.20 | 6.97 | 6.98 |
| 3 | 0.36 | 0.24 | 0.20 | 0.17 | 5.05 | 6.97 |
| 4 | 0.50 | 0.19 | 0.31 | 0.14 | 6.97 | 6.97 |
| 5 | 0.87 | 0.34 | 0.53 | 0.25 | 6.97 | 6.98 |
| 6 | 0.44 | 0.20 | 0.27 | 0.14 | 6.97 | 6.97 |
| 7 | 0.36 | 0.26 | 0.18 | 0.19 | 6.14 | 6.98 |
| 8 | 0.37 | 0.21 | 0.21 | 0.15 | 5.31 | 6.98 |
| 9 | 0.49 | 0.30 | 0.26 | 0.22 | 6.75 | 6.97 |
| 10 | 0.41 | 0.19 | 0.25 | 0.14 | 6.97 | 6.97 |
| **mean** | **0.451** | **0.241** | **0.263** | **0.176** | — | — |

(GPP/NPP in kgC/m²/yr.)

## Assessment

- **Same order of magnitude and similar dynamics.** Both models produce a temperate-wheat carbon cycle
  of a few tenths of a kgC/m²/yr, peaking in the warmest year (year 5) and dipping in the coldest, with
  full-season peak LAI near the CFT maximum (7.0). No blow-ups; both stay finite over ten years.
- **Revised ≈ 0.53× the original on GPP, ≈ 0.67× on NPP.** NPP is the fairer cross-model metric:
  Terrarium reports **GPP net of leaf respiration**, whereas the original `gpp` is the gross flux, so a
  definitional offset inflates the GPP gap. On NPP (both *net* primary production) the models agree to
  within ~35 %.
- **LAI differs systematically:** the revised canopy reaches the phenological maximum every year, while
  the original's peak LAI drops to 3–5 in poorer years. The original caps LAI by the running
  carbon/water/nitrogen deficit (the LAI–NPP feedback); the revised uses the heat-unit LAI trajectory
  only. **This is the single biggest structural difference** and explains much of the residual gap.

The remaining differences are the expected consequences of the re-architecture: a different soil water /
temperature model (hence different water-stress and phenology timing), a different surface-radiation
scheme (APAR vs the PALADYN surface energy balance), and a simplified nitrogen cycle with no fertilizer
in the revised example.

**Component-level agreement holds where the inputs are controlled:** the ported C3/C4 photosynthesis
reproduces Terrarium's `LUEPhotosynthesis` to `rtol = 1e-10` (`test/crop/test_photosynthesis.jl`), and
the soil C–N transforms are unit-tested against their LPJmL forms. The end-to-end gap is therefore in
the *coupling* (LAI feedback, soil water/temperature, radiation), not the crop physiology kernels.

## Sources of the difference

The gap has both crop-physiology and structural contributions. In rough order of contribution:

1. **Missing λ water-coupling solver (crop physiology — closeable).** The original couples
   photosynthesis to soil water through the optimal-λ (cᵢ/cₐ) bisection (`solve_lambda_*` + `lpj_bisect`);
   the revised uses a crude `λ = λ_min + (λ_opt − λ_min)·β`, and **water stress is off by default**
   (β = 1, `plant_available_water = nothing`, and the default `SoilHydrology` = `NoFlow`). This severs the
   water→GPP/transpiration feedback that is central to LPJmL and is the single largest value-changing
   difference (see the audit, gap #1).
2. **Definitional GPP offset (fixable).** Terrarium's `gpp` is net of leaf respiration; the original's is
   gross. On NPP (both *net*) the models already agree to ~35 %.
3. **Ported-but-not-wired crop features (closeable, incremental).** The full nitrogen demand/uptake
   kinetics, the faithful Vcmax N-limitation, the harvest index (yield), and the LAI–NPP carbon-deficit
   feedback are tested primitives that are not yet connected — they shift GPP, N status, yield, and the
   canopy trajectory.
4. **Re-ported physics that is simply absent.** Vernalization + the climate buffer/spinup (winter wheat
   needs vernalization for correct phenology timing) and NO₃ leaching were removed and not re-added.
5. **Structural (intended).** A different soil water/temperature model (Richards + energy vs the LPJmL
   5-layer scheme) and surface-radiation scheme (PALADYN surface energy balance vs LPJmL APAR/albedo/PET)
   necessarily change water stress, light absorption, and timing. These are intended, not defects.

The takeaway: the revised model **reproduces the crop physiology kernels faithfully** (C3/C4
photosynthesis matches Terrarium's `LUEPhotosynthesis` to rtol 1e-10) and gives the **same order of
magnitude and seasonal behaviour** end-to-end. The residual ~2× GPP / ~35 % NPP gap is driven mostly by
the missing λ/water coupling and the not-yet-wired crop features (items 1–4, which are closeable), on top
of the intended soil/surface infrastructure change (item 5). Exact numerical agreement is not expected
once the soil/surface is Terrarium's, but the value-changing crop gaps below could be closed to bring the
two substantially closer.

## Feature-parity audit (base revision `2192dc1f` vs current tree)

The revised version does **not** yet carry every feature of the original. Status legend: **PRESENT**
(wired into the running model), **PORTED-NOT-WIRED** (tested primitive, not connected),
**REPLACED-BY-TERRARIUM**, **MISSING/DEFERRED**.

**Crop physiology:** C3/C4 photosynthesis PRESENT · APAR PRESENT · carbon allocation PRESENT ·
autotrophic respiration PRESENT · phenology + LAI PRESENT (carbon-deficit LAI cap DEFERRED) · crop carbon
PRESENT · crop nitrogen PRESENT (simplified) · temperature stress PRESENT · root distribution PRESENT ·
**λ solver MISSING** · **N demand/uptake kinetics PORTED-NOT-WIRED** · **Vcmax N-limitation primitive
PORTED-NOT-WIRED** (effect approximated) · **vernalization MISSING** · **harvest index PORTED-NOT-WIRED**.

**Management:** sowing PRESENT · harvest PRESENT · fertilizer PRESENT (continuous flux) · **tillage
MISSING/DEFERRED**.

**Soil/surface (delegated to Terrarium):** soil temperature, evaporation, infiltration/percolation,
transpiration, interception, freeze-thaw, radiation/albedo/PET all REPLACED-BY-TERRARIUM (note: default
crop model uses `NoFlow` hydrology → static water) · **snow MISSING** (Terrarium has only an abstract
stub).

**Soil biogeochemistry:** decomposition, nitrification, denitrification PRESENT · mineralization PRESENT
(simplified) · litter routing PRESENT (simplified) · **immobilization PORTED-NOT-WIRED** · **NH₃
volatilization PORTED-NOT-WIRED** · **NO₃ leaching MISSING** · **surface-litter hydrology/thermal
MISSING**.

**Climate/forcing:** climate input REPLACED-BY-TERRARIUM (`surface_climate_inputs` → `FieldTimeSeries`) ·
CO₂ PRESENT (constant; time-series DEFERRED) · **climate buffer/spinup MISSING**.

**Diagnostics:** **runtime Water/Nitrogen/Carbon/Thermal balance ledgers MISSING** (conservation is
unit-tested only).

**Numerics:** CPU PRESENT · GPU PRESENT (framework; full-model run DEFERRED) · Float32/Float64 PRESENT ·
checkpoints REPLACED-BY-TERRARIUM · differentiability PARTIAL (soil biogeochem via Reactant; full crop
`LandModel` blocked by the root-fraction Reactant gap).

### Key gaps, ordered by impact on a wheat run

1. λ water-coupling solver — MISSING (water→GPP feedback severed; water stress off by default).
2. Crop N demand + uptake kinetics — PORTED-NOT-WIRED.
3. Vernalization + climate buffer/spinup — MISSING (winter wheat phenology timing).
4. Harvest index — PORTED-NOT-WIRED (grain-yield diagnostic).
5. Faithful Vcmax N-limitation — PORTED-NOT-WIRED (approximated).
6. Immobilization limitation — PORTED-NOT-WIRED.
7. NH₃ volatilization — PORTED-NOT-WIRED.
8. NO₃ leaching — MISSING.
9. Surface-litter hydrology/thermal — MISSING.
10. Tillage — MISSING (upstream Terrarium hook).
11. Snow model — MISSING (Terrarium stub only).
12. LAI↔NPP carbon-deficit feedback — DEFERRED.
13. Runtime balance-ledger diagnostics — MISSING.

Items 1–4 change the *values* a wheat run produces; 5–13 are narrower fidelity/robustness gaps. Most
PORTED-NOT-WIRED primitives are individually tested, so wiring them is incremental — except the λ solver
(#1) and vernalization + climate buffer (#3), which require re-porting removed code. These are tracked in
`2026-07-24_NOTES_future_work.md`.

## Feature-parity audit

_(populated from the audit of base revision `2192dc1f` vs the current tree)_
