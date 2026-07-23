# Agrocosm.jl — Agent Rules

Agrocosm.jl is a **downstream package built on [Terrarium.jl](../Terrarium-copy/Terrarium.jl)**
(the `[sources]` dev dependency currently points at the `Terrarium-copy` checkout on branch
`mg/adjust-for-neuralcrop`, which carries the vegetation-dispatch widenings the crop processes
depend on). It contributes
crop-specific processes (C3/C4 photosynthesis, carbon/nitrogen allocation, phenology, coupled soil
C–N biogeochemistry, and crop management) plus a managed-crop land model. Terrarium supplies **all**
infrastructure and the physical soil/surface processes.

> **Status:** mid-migration from a standalone discrete-time LPJmL-derived model onto Terrarium.
> The active plan lives in `docs/dev/2026-07/2026-07-23_PLAN_terrarium_migration.md`. Read it before
> making structural changes.

## Framework rules come from Terrarium

Terrarium's `AGENTS.md` is authoritative for framework conventions and is **required reading**:
`../Terrarium.jl/AGENTS.md`. In particular, follow its rules on:

- **Continuous-time dynamics only.** All process implementations must be well-formed ODEs/PDEs.
  Discrete-time updates are prohibited *except* clearly documented and justified exceptions. In
  Agrocosm the sanctioned exceptions are genuine discrete management events (sowing, harvest,
  tillage, residue transfer), implemented as documented Oceananigans callbacks. Everything that
  admits a continuous formulation (e.g. fertilizer as a time-distributed input flux) must be
  expressed as tendencies/inputs.
- **Kernels & GPU compatibility** (`@kernel`/`@index`, type-stable, allocation-free, no reachable
  throw paths, `ifelse` over `if`/`else`, no hardcoded `Float64`).
- **Differentiability (Enzyme) and Reactant compatibility.**
- **The process/model interface**: `variables`, `initialize!`, `compute_auxiliary!`,
  `compute_tendencies!`, optional `closure!`/`invclosure!`, implemented across the three levels
  (interface methods → `@kernel` entry points → scalar `compute_*` primitives).

## Package-specific conventions

- **Imports.** Source code uses explicit imports (`using Terrarium` for the framework surface, plus
  explicit `using SpeedyWeatherInternals.ParameterEditing: @parameterized, @param, ...` for the
  parameter macros, and explicit kernel imports). Examples/tests use `using Agrocosm` / `using
  Terrarium` and never explicitly import exported names. ExplicitImports.jl is a test dependency.
- **Parameters.** Crop and CFT parameters are defined with ModelParameters via
  `@parameterized`/`@param` (units + bounds), exposing them to `ParameterEditing.reconstruct` for
  differentiable calibration — never `Parameters.@with_kw`.
- **Docstrings.** DocStringExtensions with `$TYPEDSIGNATURES`; always `jldoctest`, never plain
  `julia` blocks.
- **Formatting.** Runic (`fredrikekre/runic-action` in CI; config-free).
- **Implementation plans.** Major changes are prefaced by a dated plan under `docs/dev/YYYY-MM/`
  following the template in Terrarium's `AGENTS.md`.

## Testing

- Test-only dependencies live in `test/Project.toml` (not the root project). Run individual test
  files by activating the test env from the project env with
  [TestEnv.jl](https://github.com/JuliaTesting/TestEnv.jl):
  ```julia
  # started with: julia --project=.
  using TestEnv; TestEnv.activate()
  include("test/processes/crop/test_photosynthesis.jl")
  ```
- Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- New crop processes require CPU + GPU tests and Enzyme adjoint tests, plus conservation checks in
  place of the legacy bespoke balance ledgers.

## Migration guardrails

- The migration is a single long-lived feature branch (`mg/revise-terrarium`); there is **no working
  end-to-end model until the crop physics are ported** (plan Phase 6). Do not expect `using
  Agrocosm` to reproduce the old simulation API mid-migration.
- When deleting legacy infrastructure, preserve embedded physics flagged in the plan (temperature
  stress, vernalization, root distribution, nitrate advective transport, the daily process ordering,
  and the flux/stock taxonomy behind the balance ledgers). Git is the journal for anything removed.
