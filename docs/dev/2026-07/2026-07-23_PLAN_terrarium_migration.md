# Re-architect Agrocosm.jl onto the Terrarium.jl framework

> Status: **in progress**. Phases 0 (framework toolchain + conventions) and 1 (infrastructure swap)
> are complete: Agrocosm is now a downstream Terrarium package whose skeleton loads on the framework,
> with the crop/soil physics retained on disk for porting in Phases 2–6.

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
> - **Physical soil/surface/climate process files are retained on disk (excluded from the module)**
>   rather than deleted in Phase 1. Every physics file depends on the deleted infrastructure
>   (`launch_1D!`, `@unpack`, legacy state containers) and cannot compile yet; keeping them in place
>   preserves the reference implementation for porting. They are deleted incrementally as each is
>   superseded by a Terrarium config/process in Phases 2–3. `src/Agrocosm.jl` carries a documented
>   manifest of these pending files.
> - **The parameter migration to `@parameterized`/`@param` is deferred.** `default_params.jl` and
>   `pft.jl` are already infrastructure-free `@kwdef` structs, so Phase 1 keeps them as-is (retaining
>   the LPJmL defaults as physics reference). Their conversion to ModelParameters calibration hooks
>   happens with the consuming processes (Phase 3) and the CFT presets (Phase 5), where it can be
>   validated against behaviour.

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
  Richards hydrology + stratigraphy) and the surface stack; delete Agrocosm's physical soil/surface
  code. Contribute upstream only genuinely-missing pieces (e.g. multi-layer snow, already on
  Terrarium's roadmap).
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
  Terrarium idioms; final dead-code and redundancy audit.

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
