# Re-architect Agrocosm.jl onto the Terrarium.jl framework

> Status: **in progress**. Phases 0–2 complete; Phase 3 (crop physiology) underway. Ported and
> unit-validated: crop root distribution, C3/C4 photosynthesis (matches Terrarium's LUEPhotosynthesis
> to rtol 1e-10), crop stomatal conductance, and crop carbon dynamics. **A first coupled crop model
> runs end-to-end on CPU and produces positive GPP** (crop photosynthesis + stomatal conductance +
> carbon dynamics assembled into a Terrarium `LandModel`, enabled by widening Terrarium's vegetation
> dispatches to the abstract interfaces). Remaining Phase 3: crop phenology (heat units), crop
> nitrogen, plant-available-water, and soil C–N biogeochemistry. Legacy physics retained on disk until
> the Phase 6 cleanup.

Date of initial draft: 2026-07-23

Base revision: 2192dc1f6c11c4283511104c2fa64860a1224bbe

## Originating prompt

> Make a plan to revise the complete repository and package to actually use Terrarium.jl as a
> framework. Inspect Terrarium.jl's AGENTS.md and documentation first. Try to reduce code
> redundancies to an absolute minimum.

Confirmed design decisions:

1. **Repository strategy — downstream package.** Agrocosm keeps its own repository and UUID,
   depends on Terrarium, uses `using Terrarium`, and adds only crop-specific processes plus a
   managed-crop land model. Terrarium supplies all infrastructure and physical soil/surface
   processes.
2. **Management actions — hybrid.** Actions that admit a continuous formulation (e.g. fertilizer
   as a time-distributed input flux) are expressed as tendencies/inputs; genuine discrete events
   (sowing, harvest, tillage) are implemented as documented Oceananigans callbacks, invoking
   Terrarium's continuous-time exception clause.
3. **Sequencing — full rewrite in one branch.** The crop-soil model is rebuilt on Terrarium in a
   single long-lived feature branch before end-to-end runs are restored.

## Revision log

> 2026-07-23: initial draft.
>
> 2026-07-23: executed Phases 0 and 1.
>
> - **Phase 0.** Reworked `Project.toml` to a downstream Terrarium package (Terrarium as a dev path
>   dependency via `[sources]`; DocStringExtensions, SpeedyWeatherInternals, Unitful; dropped the
>   standalone infra deps). Added `test/Project.toml` (TestEnv, ExplicitImports), `AGENTS.md`
>   (deferring framework rules to `Terrarium/AGENTS.md`), `.claude/CLAUDE.md`, and a Runic CI
>   workflow. Validated the toolchain with a CPU `SoilModel` column spike
>   (`docs/dev/2026-07/spike_soil_column_cpu.jl`), which runs a 3-day soil integration end-to-end.
>   The **GPU spike is deferred** — the dev host is Apple Silicon with no CUDA; GPU validation moves
>   to CI/a CUDA host.
> - **Phase 1.** Deleted the pure-infrastructure files outright (recoverable via git):
>   `utils/{kernel_launch,conversions,load_nc,visualization}.jl`, `input_output/**`,
>   `diagnostics/**`, `simulations/**`. Rewrote `src/Agrocosm.jl` to `using Terrarium`, including only
>   the infrastructure-free files (`parameters/{default_params,pft}.jl`, `numerics/lpj_bisect.jl`) and
>   re-exporting the crop parameter/CFT-registry/bisection API. `using Agrocosm` now precompiles and
>   loads; the surviving numeric + CFT-registry tests pass (19/19).
>
> Deviations from the original phase text, made for safety and to ease later porting:
>
> - **All superseded legacy source files are retained on disk (excluded from the module) until a
>   single final cleanup in Phase 6**, rather than deleted incrementally (see the "Legacy-file
>   retention policy" section). Every physics file depends on the deleted infrastructure
>   (`launch_1D!`, `@unpack`, legacy state containers) and cannot compile yet; keeping them in place
>   preserves the reference implementation for porting. `src/Agrocosm.jl` carries a documented
>   manifest of these pending files. Only the pure-infrastructure files with nothing to port were
>   removed in Phase 1.
> - **The parameter migration to `@parameterized`/`@param` is deferred.** `default_params.jl` and
>   `pft.jl` are already infrastructure-free `@kwdef` structs, so Phase 1 keeps them as-is (retaining
>   the LPJmL defaults as physics reference). Their conversion to ModelParameters calibration hooks
>   happens with the consuming processes (Phase 3) and the CFT presets (Phase 5), where it can be
>   validated against behaviour.
>
> 2026-07-23: executed Phase 2 (reuse Terrarium soil & surface).
>
> - Documented the legacy→Terrarium soil/surface/climate configuration mapping and the crop coupling
>   seams into Terrarium's vegetation stack (`AbstractPhotosynthesis`, stomatal conductance,
>   respiration, phenology, carbon dynamics, root distribution, plant-available-water) and the soil
>   `biogeochem` slot, in `docs/dev/2026-07/2026-07-23_PHASE2_soil_surface_config.md`.
> - Validated the coupled Terrarium `LandModel` (soil energy+Richards hydrology + surface energy
>   balance + surface hydrology + prescribed atmosphere + vegetation carbon) on CPU with a
>   single-timestep spike (`docs/dev/2026-07/spike_land_column_cpu.jl`): the full coupled compute path
>   runs and exposes 64 state variables (soil temperature/moisture, latent/ground heat flux, GPP, LAI,
>   …). Per the retention policy, no legacy soil/surface files were deleted.
> - Identified genuinely-missing pieces for later upstreaming/implementation: a snow scheme (no snow
>   component in Terrarium's coupled path yet), nitrate advective transport (port with soil C–N),
>   and crop canopy PAR / LAI-dependent albedo (port with crop photosynthesis).
> - Known numerical note: a sustained VanGenuchten Richards integration is stiff at Δt = 60 s under a
>   near-saturated column (the saturation prognostic diverges); selecting a stable timestepper/Δt for
>   the crop model is a Phase 5 concern. Multi-day *soil* integration is separately green
>   (`spike_soil_column_cpu.jl`, 3-day `SoilModel` run).

> 2026-07-23: started Phase 3 (port crop physiology).
>
> - Established the downstream crop-process pattern and ported the first process:
>   `CropRootDistribution <: Terrarium.AbstractRootDistribution`
>   (`src/crop/root_distribution.jl`), re-expressing LPJmL's cumulative root profile
>   `Y(d) = 1 - β^d` as a continuous, column-normalized root density, with a CPU test. Added the
>   crop-process dependencies (KernelAbstractions, Oceananigans) and module imports.
> - Sequencing finding: after the self-contained root distribution, the remaining crop physiology
>   (phenology, C3/C4 photosynthesis + λ solver, autotrophic respiration, carbon allocation, LAI,
>   crop nitrogen) and soil C–N biogeochemistry are **mutually coupled through shared crop state**
>   (~2,300+ source lines). Faithful continuous-time reformulation of these is only meaningfully
>   testable once assembled into a crop `VegetationModel`/`LandModel`. Recommended approach for the
>   remainder of Phase 3: port the coupled physiology as a cohesive unit alongside a minimal crop
>   vegetation model (bringing part of Phase 5 forward), so each process's continuous tendencies can
>   be validated end-to-end and against the legacy daily-loop numerics, rather than as isolated
>   stubs. The discrete sowing/harvest/senescence lifecycle events remain Phase 4.
> - `temp_stress` is a scalar primitive consumed by photosynthesis (not a standalone process); it
>   ports as a Level-III primitive inside the crop photosynthesis process.
>
> 2026-07-23: ported C3/C4 crop photosynthesis (plan Phase 3c).
>
> - `CropPhotosynthesis <: Terrarium.AbstractPhotosynthesis` (`src/crop/photosynthesis.jl`),
>   continuous-time C3/C4 with the full three-level structure. C3 reproduces Terrarium's
>   `LUEPhotosynthesis` scalar `(Rd, An)` to `rtol = 1e-10` across 72 input combinations
>   (`test/crop/test_photosynthesis.jl`); C4 uses `c₂ = 1`, `φ = min(1, λ/λ_mc4)`, tested for gating,
>   light response, φ saturation, and the 55/45 °C cutoffs. Design + finding in
>   `2026-07-23_PHASE3_photosynthesis_design.md`.
> - **Integration finding:** `MedlynStomatalConductance` dispatches on `LUEPhotosynthesis`
>   specifically (not `AbstractPhotosynthesis`), which blocks injecting `CropPhotosynthesis` into
>   `VegetationCarbon`. Fix is a one-line upstream Terrarium widening, or the planned crop
>   stomatal-conductance (LPJmL λ solver) port; until then the crop photosynthesis is validated at the
>   process/unit level and end-to-end assembly is deferred to that next port.
>
> 2026-07-23: ported crop stomatal conductance and assembled the first coupled crop model.
>
> - `CropStomatalConductance <: Terrarium.AbstractStomatalConductance` (`src/crop/stomatal_conductance.jl`)
>   supplies `leaf_to_air_co2_ratio` (λ = LAMBDA_OPT reduced toward a stressed floor by β) and
>   `canopy_water_conductance`; dispatches on `AbstractPhotosynthesis`, so it pairs with
>   `CropPhotosynthesis` (avoiding the Terrarium Medlyn edit). Unit-tested
>   (`test/crop/test_stomatal_conductance.jl`). The full LPJmL supply=demand λ bisection remains a
>   documented refinement.
> - **End-to-end assembly validated:** `spike_crop_vegetation_model.jl` injects both crop processes
>   into `VegetationCarbon`/`LandModel` and takes a coupled CPU timestep — the crop physiology runs in
>   the full stack (λ = 0.8 well-watered). GPP is still 0 because the PALADYN-default phenology/carbon
>   slot feeds a non-physical `leaf_area_index` (≤ 0); meaningful GPP awaits the crop LAI
>   (phenology/carbon, P3d) and crop PAW (P3f) ports. Added defensive guards so the crop processes are
>   robust to the out-of-range upstream values from the not-yet-ported slots: β clamped to [0,1] in
>   both crop kernels, and canopy cover clamped to ≥ 0 (a negative default LAI was otherwise driving
>   `canopy_water_conductance` negative).
>
> 2026-07-23: ported crop carbon dynamics; key architectural finding for Phase 5.
>
> - `CropCarbonDynamics <: Terrarium.AbstractVegetationCarbonDynamics` (`src/crop/carbon_dynamics.jl`):
>   prognostic `carbon_vegetation`, `balanced_leaf_area_index = C_veg/(2/SLA + awl)` (clamped ≥ 0), and
>   a **dimensionally-correct per-second** turnover tendency (Terrarium's `PALADYNCarbonDynamics`
>   carries turnover as yr⁻¹ but the timestepper integrates per-second — its own flagged TODO — which
>   collapses `carbon_vegetation` negative in a single step, the true source of the negative default
>   LAI). Unit-tested (`test/crop/test_carbon_dynamics.jl`).
> 2026-07-23: completed crop nitrogen cycle (P3e) and advanced soil C–N (P3f).
>
> - **P3e complete.** Crop nitrogen allocation (`CropNitrogenAllocation` / `allocate_crop_nitrogen`,
>   `src/crop/nitrogen_allocation.jl`, LPJmL `crop_nitrogen`): total plant N redistributed among
>   leaf/root/storage/pool by carbon-to-target-C:N weights, conserving N. Unit-tested. Together with
>   the earlier demand, NO₃/NH₄ uptake kinetics, and Vcmax limitation, the crop nitrogen cycle is
>   ported as tested scalar physics.
> - **P3f soil C–N.** `CropSoilCarbon` (`src/crop/soil_carbon.jl`, LPJmL `soil_carbon`): first-order
>   pool decomposition `(1 − exp(−rate·response))·pool`, litter routing to fast/slow/atmosphere
>   (carbon-conserving), and heterotrophic respiration. `CropNitrification`
>   (`src/crop/nitrification.jl`, LPJmL `nitrogen_transform`): gross NH₄→NO₃ nitrification as a peaked
>   WFPS moisture factor × Gaussian temperature factor × atan pH factor, capped at the ammonium stock,
>   with the N₂O split. `CropDenitrification` (`src/crop/denitrification.jl`, LPJmL
>   `nitrogen_transform`): gross NO₃ denitrification as a WFPS moisture factor × an organic-carbon
>   availability factor (fast + slow carbon with a peaked soil-temperature response) × nitrate, capped
>   at the nitrate stock, split into N₂O/N₂. All unit-tested (incl. the nitrification moisture and
>   denitrification temperature/moisture closed forms). Note: the nitrification/denitrification
>   moisture/temperature factors use non-negative bases before their non-integer powers so they are
>   throw-free (kernel/Reactant-safe), matching the LPJmL "0 outside support" behaviour.
> - Completed the soil-N transform set: `CropVolatilization` (`src/crop/volatilization.jl`, NH₃ flux =
>   aqueous-NH₃ fraction × Henry's constant × wind mass-transfer, capped at top-layer ammonium) and
>   `CropNitrogenMineralization` (`src/crop/mineralization.jl`, the `litter_C/soil_CN − litter_N`
>   immobilization demand and its Michaelis-Menten limitation by available mineral N). Both
>   unit-tested. **P3f soil C–N is now ported as tested scalar physics** (carbon pool decomposition +
>   litter routing, nitrification, denitrification, volatilization, mineralization/immobilization);
>   the multi-layer pool state, C-shift/litter routing between pools, and the coupled assembly to the
>   soil hydrology/temperature are Phase 5.
>
> 2026-07-23: completed crop carbon physics (maintenance respiration + organ allocation).
>
> - `CropMaintenanceRespiration` (`src/crop/maintenance_respiration.jl`, LPJmL crop `respiration`):
>   per-organ maintenance respiration `carbon·respcoeff·k·nc_ratio·f(T)` with a Lloyd-Taylor
>   temperature response (root at soil temperature, storage/pool at air temperature). With growth
>   respiration this completes crop autotrophic respiration `Ra = Rm + Rg`.
> - `CropCarbonAllocation` (`src/crop/carbon_allocation.jl`, LPJmL/SWAT `carbon_allocation`): the
>   SWAT-style root fraction of biomass (declines through the season, rises under water/N stress) and
>   the LAI/SLA-constrained leaf carbon. Both unit-tested. The crop carbon-cycle physics
>   (photosynthesis → respiration → NPP → allocation → LAI, with the harvest index) is now ported as
>   tested scalar primitives.
>
> 2026-07-23: ported crop harvest index + growth respiration (carbon allocation, part of P3d).
>
> - `CropHarvestIndex` / `crop_harvest_index` (`src/crop/harvest_index.jl`, LPJmL `carbon_allocation`):
>   the storage-organ carbon fraction — a phenology-driven optimum (sigmoid in `fphu`) scaled between
>   `himin` and `hiopt`, reduced by the water-deficit sufficiency factor, with HI > 1 scaled about 1
>   for above-ground:total-biomass crops. Unit-tested against the closed form, water-stress limits,
>   phenology dependence, and the HI > 1 case (`test/crop/test_harvest_index.jl`).
> - `CropGrowthRespiration` (`src/crop/growth_respiration.jl`): the LPJmL `r_growth = 0.25`
>   post-maintenance carbon split — `Rg = r_growth·max(0, GPP − Rm)`, `NPP = GPP − Rm − Rg`.
>   Unit-tested.
>
> 2026-07-23: ported crop plant-available-water stress (start of P3f).
>
> - `soil_moisture_limiting_factor` and `plant_available_water` (`src/crop/plant_available_water.jl`):
>   the water-stress factor β = clamp((θ − θ_wilting)/(θ_field_capacity − θ_wilting), 0, 1) that crop
>   photosynthesis/stomatal conductance consume (currently defaulted to 1), plus per-layer
>   plant-available water. Unit-tested (`test/crop/test_plant_available_water.jl`). The depth-integrated,
>   root-weighted coupling to Terrarium's soil hydraulics is wired in the crop vegetation model
>   (Phase 5); soil C–N biogeochemistry is the rest of P3f.
>
> 2026-07-23: ported crop leaf-nitrogen limitation of Vcmax (part of P3e).
>
> - `CropNitrogenVcmaxLimit` + `nitrogen_limited_vcmax` (`src/crop/nitrogen_limitation.jl`): LPJmL's
>   `limit_vcmax_by_nitrogen` as a tested scalar primitive — structural leaf N (`ncleaf_min·leaf_C`) is
>   protected, only the excess supports Rubisco, temperature-scaled, capping the potential Vcmax and
>   returning the retained fraction ∈ [0,1]. Unit-tested against the closed form plus the non-limiting,
>   structural-protection, zero-potential, and temperature-dependence cases
>   (`test/crop/test_nitrogen_limitation.jl`). Applied to the crop photosynthesis Vcmax within the crop
>   nitrogen coupling (Phase 5); the full crop N cycle (demand/uptake/allocation) is the rest of P3e.
>   Also ported the complementary **crop nitrogen demand** (`CropNitrogenDemand` /
>   `crop_nitrogen_demand`, `src/crop/nitrogen_demand.jl`, LPJmL `ndemand_crop`): leaf demand = the
>   Rubisco N requirement implied by Vcmax (the inverse of the Vcmax limitation — a cross-check test
>   confirms feeding leaf demand back through the limitation recovers the Vcmax) plus structural leaf
>   N, and total demand adds root/pool/storage at the leaf N:C ratio scaled by organ C:N ratios.
>   Unit-tested against the closed form, the demand↔limitation inversion, and temperature/carbon
>   dependence (`test/crop/test_nitrogen_demand.jl`).
>   Also ported the **root nitrogen uptake kinetics** (`CropNitrogenUptakeKinetics`,
>   `nitrogen_uptake_temperature_response`, `root_nitrogen_uptake_potential`,
>   `src/crop/nitrogen_uptake.jl`, LPJmL `nuptake_crop`): per-pool Michaelis-Menten uptake
>   `vmax·(kmin + N/(N + Km·scale))·root_factor` capped at the available N, with the LPJmL parabolic
>   soil-temperature response (normalized to 1 at the reference temperature). Unit-tested
>   (`test/crop/test_nitrogen_uptake.jl`). Still remaining for P3e: nitrogen allocation across organs
>   and wiring demand/uptake/limitation into the coupled model (needs soil mineral-N pools from P3f).
>
> 2026-07-23: ported crop phenology (LAI heat-unit trajectory).
>
> - `CropPhenology <: Terrarium.AbstractPhenology` (`src/crop/phenology.jl`): the LPJmL
>   heat-unit-driven leaf-area-index trajectory `LAI = f(fphu)·laimax` — a logistic-like rise to the
>   senescence onset `fphusen` then a power-law decline to a harvest floor. Unit-tested against the
>   trajectory's exact design anchors (`flaimaxc` at `fphuc`, `flaimaxk` at `fphuk`, peak 1 at
>   `fphusen`, bare at `fphu = 1`) plus monotonicity and bounds (`test/crop/test_phenology.jl`). The
>   heat-unit fraction `fphu` is supplied as the input `phenology_heat_unit_fraction`; its prognostic
>   accumulation `d(HU)/dt = max(0, T − T_base)` (with vernalization/photoperiod modifiers, from the
>   legacy `phenology`/`climbuf`) is integrated in the crop vegetation model (Phase 5). This replaces
>   the carbon-pool-equilibrium LAI with the faithful growing-season LAI shape.
>
> 2026-07-23: first crop coupled model producing positive GPP.
>
> - Resolved the sibling-coupling by **widening Terrarium's vegetation dispatches to the abstract
>   interfaces** in `../Terrarium-copy/Terrarium.jl` (branch `mg/adjust-for-neuralcrop`, commit
>   "Widen vegetation-process dispatch to abstract interfaces"): `MedlynStomatalConductance` →
>   `AbstractPhotosynthesis` and `PALADYNAutotrophicRespiration` → `AbstractVegetationCarbonDynamics`
>   (strictly more general; LUE/PALADYN still match; each reads only fields present on every
>   implementation). Agrocosm's `[sources]` now points at the Terrarium-copy checkout.
> - With those widenings, `spike_crop_vegetation_model.jl` assembles a `LandModel` whose
>   `VegetationCarbon` uses crop photosynthesis + crop stomatal conductance + crop carbon dynamics,
>   with the (coupling-free) `PALADYNPhenology` and the now-compatible `PALADYNAutotrophicRespiration`,
>   and `vegetation_dynamics = nothing` (a single crop cell needs no PFT-spreading dynamics). It runs
>   end-to-end on CPU and **produces positive GPP** (≈ 1.5e-8 kgC/m²/s ≈ 1.3 gC/m²/day),
>   `net_assimilation`/`leaf_respiration` finite and positive, λ = 0.8, canopy conductance positive.
>   This is the first coupled crop model on the Terrarium stack producing carbon assimilation.
> - `PALADYNVegetationDynamics` also calls the carbon dynamics' `compute_λ_NPP` method (not just
>   dispatches on its type), so it is not made compatible by dispatch-widening alone; it is bypassed
>   (`nothing`) for the single-cell crop rather than widened. A crop `VegetationModel` (Phase 5) can
>   wire it if PFT dynamics are needed.
>
> - **Architectural finding (drives Phase 5):** Terrarium's PALADYN vegetation processes are
>   concretely coupled to each other's *types*, not the abstract interfaces —
>   `MedlynStomatalConductance` dispatches on `LUEPhotosynthesis`, and
>   `PALADYNAutotrophicRespiration.compute_autotrophic_respiration`/`compute_Ra` dispatch on
>   `PALADYNCarbonDynamics` (`respiration/autotrophic_respiration.jl:133,185`). Swapping any single
>   crop process into `VegetationCarbon` therefore cascades into its PALADYN siblings. Conclusion: the
>   coupled crop physiology must be assembled in a **dedicated crop `VegetationModel`** (plan Phase 5)
>   that wires the crop processes together, rather than by slot-swapping into `VegetationCarbon`. The
>   `spike_crop_vegetation_model.jl` therefore swaps only the photosynthesis + stomatal-conductance
>   slots (which pair cleanly); `CropCarbonDynamics` is delivered and unit-validated, ready for that
>   crop vegetation model. (An alternative is to widen the offending Terrarium dispatches to the
>   abstract supertypes upstream — a small, correct framework improvement.)

> 2026-07-23: started Phase 5 (crop vegetation model assembly).
>
> - `CropPhenologyDynamics` (`src/crop/phenology_dynamics.jl`): the prognostic phenological-heat-unit
>   accumulator — `d(HU)/dt = max(0, T_air − T_base)/seconds_per_day` (continuous-time replacement for
>   LPJmL's daily heat-unit sum), producing the heat-unit fraction `fphu` that drives the LAI
>   trajectory. Unit-tested.
> - `CropVegetation` (`src/crop/vegetation.jl`): a minimal crop `AbstractVegetation` wiring
>   `phenology_dynamics → phenology → stomatal_conductance → photosynthesis` with the correct
>   auxiliary/tendency ordering, slotting into a Terrarium `LandModel`. The spike
>   `spike_crop_vegetation_phenology.jl` validates it end-to-end on CPU: **heat units accumulate over
>   the integration (0 → 1.74 °C·d in 20 steps at 25 °C), and a mid-season heat-unit state
>   (fphu = 0.5) yields LAI ≈ 6.8 and positive GPP (~8 gC/m²/day)** — the crop LAI now responds to the
>   growing season (accumulated heat units) rather than a carbon-pool equilibrium.
> - **CFT presets** (`src/crop/cft_presets.jl`): `CropVegetation(NF, crop_pft(:maize))` and the
>   per-component `CropPhenology/CropPhotosynthesis/CropPhenologyDynamics(NF, pft)` constructors map the
>   12-CFT trait registry (Phase 1's `pft.jl`) onto the crop processes — the LAI trajectory, C3/C4
>   pathway + temperature thresholds, and base temperature. Unit-tested (wheat→C3, maize→C4). The crop
>   model is now configurable per crop type.
> - **Root distribution + optional plant-available-water** added to `CropVegetation`. Default is a
>   well-watered crop (β=1, robust); passing `plant_available_water=FieldCapacityLimitedPAW(NF)` couples
>   β to soil water, but **requires a clay-bearing soil texture** — the default pure-sand texture makes
>   field_capacity==wilting_point so β is NaN (documented finding).
> - **Prognostic carbon pool** (`src/crop/carbon.jl`): `CropCarbon` closes the crop carbon loop —
>   biomass (kgC/m²) accumulates NPP and partitions into leaf/root/storage each step via the ported
>   allocation + respiration primitives (`root = root_allocation_fraction(fphu, df)·biomass`,
>   `leaf = min(LAI/SLA, biomass − root)`, `NPP = GPP − Rm − Rg`, `d(biomass)/dt = NPP`), all in
>   kgC/m²/s. Wired into `CropVegetation` and the CFT presets (SLA from the registry). Spike shows
>   GPP → NPP (~75 % of GPP) and biomass accumulating over the integration; unit-tested (organ
>   conservation, NPP < GPP, water-stress root shift). The crop model now carries **two prognostics
>   (heat units + biomass)** and produces a growing-season carbon budget.
> - Remaining Phase 5: multi-pool nitrogen (leaf/root/storage) coupled to the crop N demand/uptake and
>   the soil mineral-N state; couple the soil C–N transforms to the soil biogeochemistry slot; the
>   LAI-feedback carbon deficit; per-CFT heat-unit requirement (climate-derived); optional multi-crop
>   tiling.

## Problem description

Agrocosm.jl (~10,165 source lines) is a self-contained, LPJmL-derived crop-soil model with a
**discrete daily time loop** (`daily_crop_C3!`/`daily_crop_C4!`), bespoke state containers
(`ModelState`, `Crop`, `Soil` with hand-written lifecycle selectors in
`src/simulations/model_runtime.jl`), manual index-vector domains with fixed five-layer soil
arrays, a custom `device`/kernel-launch/CUDA/Adapt backend layer, `Parameters.@with_kw`
parameters, custom NetCDF/JLD2 I/O, custom balance diagnostics, and a custom checkpoint format.

Terrarium.jl is a continuous-time, Oceananigans-based land-modeling framework that already
provides every one of those infrastructure layers plus differentiable (Enzyme) and Reactant
compilation. Terrarium mandates continuous-time ODE/PDE dynamics and prohibits discrete-time
updates except as clearly documented exceptions (`Terrarium/AGENTS.md`).

Consequently "use Terrarium as a framework" is a re-architecture, not a mechanical refactor. The
redundancy reduction is achieved by deleting Agrocosm's infrastructure and physical soil/surface
code (roughly half the codebase) and re-expressing the crop-specific physics on Terrarium's
`AbstractProcess`/`AbstractModel` interface.

## Background

### Terrarium interfaces the migration targets

- **Grids & Fields.** `ColumnGrid`/`ColumnRingGrid` over an `AbstractArchitecture`
  (`CPU()`/`GPU()`/`ReactantState()`); Oceananigans `Field`s with `XY()`/`XYZ()` dims.
- **State.** `variables(::Process)` returns symbolic `prognostic`/`auxiliary`/`input` variables;
  `StateVariables` allocates `Field`s, auto-creates tendency fields, and merges/promotes
  duplicates (see `docs/src/extending/state_variables.md`).
- **Process/model interface.** `variables`, `initialize!`, `compute_auxiliary!`,
  `compute_tendencies!`, and optional `closure!`/`invclosure!`. Three levels: top-level interface
  methods → `@kernel` entry points → scalar `compute_*(i,j[,k],grid,fields,::Process,args...)`
  primitives (see `docs/src/extending/implementing_processes.md`).
- **Models.** `AbstractModel{NF,Grid}`; hierarchy incl. `AbstractSoilModel`,
  `AbstractVegetationModel`, `AbstractLandModel`. `grid` positional; `@parameterized @kwdef`
  structs with `@component` fields; coupled models forward each interface call to their
  processes (`src/models/coupled/land_model.jl`).
- **Runtime.** `initialize(model; inputs, boundary_conditions, ...) → ModelIntegrator`;
  `run!(integrator; period=Day(10))` / `timestep!`; timesteppers `ForwardEuler`, `Heun`, `IMEX`.
- **Parameters.** `ModelParameters.@parameterized`/`@param` with `units`/`bounds`;
  `ParameterEditing.reconstruct(model, params)` gives differentiable parameter editing — directly
  serving Agrocosm's gradient-based-calibration goal.
- **I/O.** `InputSources`/`InputProvider`; `TerrariumRastersExt` for gridded inputs; Oceananigans
  output writers.
- **Extensions.** `TerrariumCheckpointingExt` (also AD reverse-pass checkpointing via
  `Reactant.Periodic`), `TerrariumReactantExt`.
- **Existing physical processes to reuse wholesale.** Soil energy (`SoilThermodynamics` with
  enthalpy/freeze closure), soil hydrology (`SoilHydrology(RichardsEq())`), soil stratigraphy /
  texture / porosity, surface energy balance, radiative fluxes, albedo, evapotranspiration
  (bare-ground & canopy), canopy interception, surface runoff, aerodynamics, prescribed
  atmosphere, physical constants, unit conversions.

### Crop-specific physics Terrarium lacks (Agrocosm must contribute as new processes)

C3/C4 mechanistic photosynthesis with the λ water-coupling bisection solver
(`src/processes/crop/lambda_solver.jl`, `src/numerics/lpj_bisect.jl`); crop carbon allocation to
organs + harvest index; crop phenology (phenological heat units / vernalization / prescribed
calendar); full crop nitrogen (demand, uptake, allocation, Vcmax limitation); coupled soil C–N
decomposition (mineralization, nitrification, denitrification, volatilization, leaching, litter
routing, surface litter); management (sowing, fertilizer, manure, tillage, harvest, residue); the
12-CFT parameter registry.

## Redundancy-elimination map

**Delete (Terrarium provides the replacement):**

| Agrocosm subsystem | Terrarium replacement |
|---|---|
| `simulations/model_runtime.jl`, `processes/initialization/**` builders, `init_states.jl` | `variables()` + `StateVariables` + `AbstractInitializer` |
| `simulations/simulation_api.jl`, `simulations/daily_crop_C3.jl`/`C4.jl` | `ModelIntegrator`, `initialize`, `run!`/`timestep!` |
| `utils/kernel_launch.jl`, `utils/conversions.jl`, CUDA/Adapt plumbing | `launch!`, `on_architecture`, architecture types |
| `parameters/default_params.jl`, `Parameters.jl` | `@parameterized`/`@param`, `PhysicalConstants`, `Unitful` |
| `input_output/**`, `utils/load_nc.jl` | `InputSources`, `TerrariumRastersExt`, Oceananigans writers |
| `processes/soil/{soil_temp,water_ice_pools,soil_water,infil_perc,pedotransfer,evaporation}.jl`, snow | `SoilThermodynamics`, `SoilHydrology(RichardsEq)`, stratigraphy/texture/porosity, `BareGroundEvaporation`, snow scheme |
| `processes/crop/{radiation,albedo,interception,transpiration}.jl` | `SurfaceEnergyBalance`, `RadiativeFluxes`, `Albedo`, `CanopyInterception`, `Evapotranspiration`, `SurfaceRunoff` |
| `processes/climate/**` | `PrescribedAtmosphere`, `Aerodynamics`, `InputSources` |
| `diagnostics/**` | `diagnostics/debugging.jl` + conservation tests |
| checkpoint code in `simulation_api.jl` | `TerrariumCheckpointingExt` |
| `utils/visualization.jl`, `Plots` dep | Oceananigans output + downstream plotting |

**Keep and port as new Terrarium `AbstractProcess` implementations:** crop C3/C4 photosynthesis +
λ solver, carbon allocation + harvest index, crop phenology, crop nitrogen, soil C–N
biogeochemistry, management events, CFT registry.

## Legacy-file retention policy

Superseded legacy source files are **retained on disk (excluded from the module) until a single
final cleanup in Phase 6**, rather than deleted in the phase that supersedes them. They are the
reference implementation used while porting physics, and keeping them in place avoids `git show`
round-trips during Phases 2–5. Only the pure-infrastructure files with nothing to port were removed
in Phase 1. Every phase below that says "supersede" leaves the old files in place; the mass deletion
of all superseded legacy source happens once, in Phase 6, after the new stack reproduces the
acceptance example.

## Summary of changes (planned phases)

Because the migration is a single-branch full rewrite, all phases land on one long-lived
feature branch; end-to-end runs are restored at Phase 6.

- **Phase 0 — Setup & conventions.** `dev` Terrarium; adopt ExplicitImports,
  DocStringExtensions (`$TYPEDSIGNATURES`), `jldoctest`, ModelParameters, split `test/Project.toml`,
  Runic, and this `docs/dev/` plan convention. Spike: run a Terrarium `SoilModel` column on
  CPU and GPU to validate the toolchain.
- **Phase 1 — Infrastructure swap.** Replace grid, state, timestepping, backend, parameters, I/O,
  constants, diagnostics, and checkpointing with Terrarium equivalents (first table block). Rewrite
  `src/Agrocosm.jl` to `using Terrarium` and re-export the crop API.
- **Phase 2 — Reuse Terrarium soil & surface.** Configure `SoilEnergyWaterCarbon` (energy +
  Richards hydrology + stratigraphy) and the surface stack for the crop context; the physical
  soil/surface code is superseded but retained per the retention policy (deleted in Phase 6).
  Contribute upstream only genuinely-missing pieces (e.g. multi-layer snow, already on Terrarium's
  roadmap).
- **Phase 3 — Port crop physiology.** New processes under Terrarium vegetation & soil-biogeochem
  abstract types: C3/C4 photosynthesis (alternative `AbstractPhotosynthesis` + λ solver),
  autotrophic respiration, carbon dynamics/allocation, LAI/canopy, root distribution, plant
  available water; crop nitrogen (demand/uptake/allocation/Vcmax limit); crop phenology; coupled
  soil C–N decomposition (new `AbstractSoilBiogeochemistry`). Each implemented level I→II→III with
  continuous-time tendencies, `variables()`, and GPU + Enzyme tests.
- **Phase 4 — Management (hybrid).** Continuous-where-feasible inputs/fluxes (e.g. fertilizer
  application distributed in time) plus documented Oceananigans callbacks for true discrete events
  (sowing, harvest, tillage, residue transfer). Each discrete exception documented and justified per
  `Terrarium/AGENTS.md`.
- **Phase 5 — Crop model & CFT registry.** A managed-crop `LandModel`/`CropModel` wiring Terrarium
  soil+surface+atmosphere with the Agrocosm crop processes; the 12 CFT parameter sets as
  ModelParameters presets; multi-crop via the planned `TiledVegetationModel`.
- **Phase 6 — Validation, AD, docs, cleanup.** Reproduce the 10-year wheat GPP/NPP example on the
  new stack; conservation, Enzyme adjoint, and Reactant global tests; rewrite README/docs/examples to
  Terrarium idioms; **delete all superseded legacy source files retained during Phases 1–5** (the
  single mass cleanup per the retention policy) and do a final dead-code and redundancy audit.

## Testing and verification

- Port the existing CPU/GPU process, simulation, and balance tests to the Terrarium test layout
  (`test/Project.toml`, `TestEnv`), replacing bespoke balance ledgers with conservation checks.
- Add Enzyme adjoint tests for each new crop process (per `test/differentiability` conventions in
  Terrarium) and a Reactant global-run smoke test.
- Acceptance: the 10-year single-cell wheat GPP/NPP trajectory reproduces the current
  `examples/` result within an agreed tolerance, and conservation residuals stay bounded.

## Documentation changes

Rewrite README, docs, and `examples/` to Terrarium idioms (`ColumnGrid`, `initialize`, `run!`).
Add process doc pages following Terrarium's Overview/Implementations/Methods/Kernel-functions
structure with `jldoctest`/`@example` blocks.

## Known limitations

- Discrete management events remain documented exceptions to the continuous-time mandate.
- Mechanistic C3/C4 photosynthesis and full crop N/soil-CN are new to the Terrarium ecosystem and
  will need review for GPU type-stability and differentiability.
- The full-rewrite sequencing means no working end-to-end model until Phase 6.

## Future work

- Upstream generally useful crop/biogeochem processes into Terrarium.
- Multi-crop stands, rotations, and dynamic sowing via `TiledVegetationModel`.
- End-to-end differentiable calibration using `ParameterEditing.reconstruct` and Reactant
  reverse-mode AD with checkpointing.
