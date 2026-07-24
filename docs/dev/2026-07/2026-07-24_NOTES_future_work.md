# Agrocosm — notes, to-dos, and future work

> Living document. Captures known limitations, deferred items, and future directions after the
> re-architecture onto Terrarium.jl. The dated migration plan
> (`2026-07-23_PLAN_terrarium_migration.md`) is the historical record; this file is where new ideas and
> open threads go.

Date started: 2026-07-24

## Known limitations / robustness

- **Reactant compatibility of the crop `LandModel` (upstream Terrarium gap).** The full crop
  `LandModel` cannot yet be built/compiled on `ReactantState`: the **static root-fraction field** fails
  to trace with `MethodError: no method matching exp(::Reactant.TracedRArray{Float32, 3})`. The
  normalization `sum(R, dims = 3)` over the `∂R∂z` `FunctionField` (whose closure calls the scalar `exp`
  in `root_density`) evaluates that closure on the whole traced z-array instead of per element.
  - This is a **Terrarium** gap, not Agrocosm-specific: Terrarium's own
    `StaticExponentialRootDistribution.root_fraction` uses the identical pattern and fails the same way;
    Agrocosm's `CropRootDistribution.crop_root_fraction` mirrors it. Terrarium's Reactant suite exercises
    only `SoilModel` configs (no vegetation), so this path is untested upstream.
  - A minimal **pure-Terrarium repro** is saved at `docs/dev/2026-07/reactant_rootfraction_repro.jl`,
    ready to drop into Terrarium's `test/reactant/` when fixing it upstream (e.g. compute the normalized
    static profile with a Reactant-traceable reduction, or materialize it on the host at init).
  - Because of this, the Reactant correctness/autodiff tests target the crop **soil biogeochemistry** (a
    `SoilModel`, no vegetation/root distribution), which compiles cleanly. *(User decision 2026-07-24:
    keep the Reactant work limited for now; this is a future to-do.)*

- **Surface energy balance stiffness (CPU runtime).** Separately, the full crop `LandModel`
  `DomainError`s **at run time on CPU** for larger timesteps: a `log` of a negative saturation in
  Terrarium's turbulent-flux thermodynamics (`canopy_evapotranspiration` → `thermodynamics`) when the
  skin temperature/humidity leaves its physical range. The managed-crop `Simulation` and the 10-year
  validation therefore run at **Δt = 600 s** with full atmospheric forcing set; larger steps overshoot
  and throw. A `clamp`/guard on the saturation argument upstream would harden this.

- **Deterministic forcing → limit cycle.** The 10-year wheat validation uses a repeating synthetic
  climate, so years 2–10 are identical. Real (noisy, trended) forcing is needed to exercise
  interannual variability.

## Deferred crop/soil features (physics ported or straightforward; wiring pending)

- **Tillage.** Modifies the topsoil bulk density / hydraulics, which Terrarium's soil stratigraphy
  owns. Needs an upstream hook (a mutable, event-updatable stratigraphy/bulk-density field) before it
  can be wired as a sowing-day management event. Legacy physics: `tillage_hydraulics!` (mixing toward a
  target bulk density by a mixing efficiency).
- **NH₃ volatilization.** The scalar physics (`ammonia_volatilization`, `CropVolatilization`) is ported
  and tested but not yet wired as a top-layer **surface** N flux out of the ammonium pool.
- **LAI ↔ NPP feedback (carbon deficit).** The fully LPJmL-faithful LAI caps the phenological LAI by a
  carbon-availability deficit; currently the LAI is the phenological trajectory only. Needs the
  running carbon-deficit state and couples naturally to sowing.
- **Per-CFT heat-unit requirement.** `heat_unit_requirement` is currently a single default;
  in LPJmL it is climate-derived (the sowing→maturity heat sum for the site), which needs the sowing
  date and a climate buffer.
- **Immobilization limitation on mineralization.** `CropNitrogenMineralization` (immobilization demand
  / limitation) is ported and tested but the soil biogeochemistry currently uses net mineralization at
  the soil C:N ratio directly; wire the immobilization limitation into the NH₄ tendency.

## Explicitly out of scope (future model generation)

- **Multi-crop tiling / rotations.** Legacy Agrocosm is single-crop; Terrarium's `TiledVegetationModel`
  is planned-but-unimplemented upstream. Omitting tiling matches the original scope. Revisit once
  Terrarium ships tiled vegetation.

## Differentiability / performance

- **Full-model AD on GPU.** Once the surface throw is fixed, differentiate the whole crop `LandModel`
  (not just the soil column) through Reactant, and benchmark on a GPU / global grid.
- **CPU Enzyme rollout AD** of the crop soil biogeochemistry (with Checkpointing.jl) compiles very
  slowly for the C–N kernels — practical AD uses the Reactant path. If a CPU rollout example is wanted,
  profile the Enzyme reverse-pass compile of the `^`-heavy nitrification/denitrification kernels.
- **Parameter calibration.** The `@parameterized` crop parameters expose bounds/units to
  `ParameterEditing.reconstruct`; build a gradient-based calibration example (cultivar/physiology
  parameters vs observed GPP/LAI/yield). See the scalar-parameter-vs-grid `NF` coupling caveat in the
  session memory before differentiating scalar physical parameters under Reactant.

## Validation / science

- Drive the model with **real gridded climate** (replacing the synthetic seasonal forcing) and compare
  GPP/NPP/yield against observations or the LPJmL reference.
- Add more CFTs to the validation (maize/C4, rice, soybean) and check the C3-vs-C4 pathway split.
- Restore **conservation-ledger** style diagnostics (carbon/nitrogen closure) as end-to-end tests, in
  the spirit of the legacy balance ledgers.

## Docs / housekeeping

- Expand the documentation beyond the getting-started + API pages: process concept pages (photosynthesis
  + λ solver, phenology, the soil C–N cycle, management), following Terrarium's doc-page conventions
  (Overview / Implementations / Methods / Kernel functions), with `@example` blocks.
- Confirm Runic formatting in CI across the new files.
- `lib/AgrocosmData/` (the data-loading companion subpackage) was left untouched by the migration;
  review whether it is still needed or should be updated to the Terrarium input path.
