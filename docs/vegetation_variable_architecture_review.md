# Vegetation variable architecture review

This note records the design review that precedes Agrocosm's
`state--flux--auxiliary--workspace` audit. It focuses on the current vegetation
variable architecture of Terrarium.jl and on keeping Agrocosm
couplable without making either framework a hard dependency.

Reviewed snapshots:

- Terrarium.jl: `4f42508179cb32fbbe37562b2a59d5036bea0a4d` (2026-07-14)
- ClimaLand.jl: `d7245431e27e52fe8f9120fe5a376e305bebf2ce` (v1.10.3, 2026-07-17)

## Terrarium.jl

Terrarium's `StateVariables` has six relevant partitions:

1. `prognostic`: fields that characterize the persistent system state;
2. `tendencies`: one field for every prognostic variable, integrated by the
   time stepper;
3. `auxiliary`: fields derivable from the current prognostic state and inputs,
   without requiring their own value at the previous time step;
4. `inputs`: spatially explicit, static or time-varying forcing fields;
5. `namespaces`: recursively nested state-variable groups;
6. `timestepper_cache`: storage owned by the numerical time stepper.

Processes declare symbolic variables with `prognostic`, `auxiliary`, and
`input`. Terrarium then allocates backend-specific Oceananigans fields and
collects declarations from coupled processes. A prognostic declaration also
creates its tendency automatically. Duplicate declarations are merged only
when dimensions and units agree; prognostic/auxiliary conflicts are rejected.

The vegetation implementation is process-composed rather than one large state
struct. `VegetationCarbon` contains photosynthesis, stomatal conductance,
respiration, phenology, carbon dynamics, vegetation dynamics, root
distribution, and plant-available-water components. In its present natural
vegetation example:

- vegetation carbon and vegetation-area fraction are prognostic;
- GPP, NPP, respiration, LAI, conductance, stress factors, and root fractions
  are auxiliary;
- NPP can be declared as an input by the carbon-dynamics process and is promoted
  to the auxiliary field produced by another component;
- `compute_auxiliary!` has an explicit dependency order;
- `compute_tendencies!` changes only prognostic variables.

Terrarium does not separately type fluxes, scientific diagnostics, and kernel
scratch arrays. All three can reside in `auxiliary`, although lazy Oceananigans
operators allow some purely derived quantities to avoid their own storage.

## ClimaLand.jl

ClimaLand presents a coarser top-level split:

- `Y`: prognostic variables advanced by the time integrator;
- `p`: auxiliary/cache variables recomputed or updated during a time step.

Its `CanopyModel` is also composed from named component models: radiation,
photosynthesis, conductance, soil-moisture stress, plant hydraulics, canopy
energy, fluorescence, respiration, and biomass. Each component declares
`prognostic_vars`, `auxiliary_vars`, their types, and their spatial domain.
The canopy aggregates these declarations into a component-nested state.

Examples illustrate the split:

- plant liquid content and prognostic canopy temperature belong to `Y`;
- water potential, root water flux, GPP, net assimilation, respiration,
  conductance, and stress factors belong to `p`;
- cache fields are deliberately preallocated because they may be expensive,
  reused, or allocation-heavy;
- boundary and inter-component fluxes are updated through a distinct
  `make_update_boundary_fluxes` stage, but are still stored in `p`;
- user-facing diagnostics are an output layer and may read either prognostic or
  cache variables. In ClimaLand terminology, "diagnostic output" is therefore
  not synonymous with "cache variable".

ClimaLand shows that a performant implementation may sometimes fuse updates
that are conceptually separate. Its documentation explicitly notes cases where
auxiliary plant-hydraulic variables are updated inside a tendency pass to avoid
another traversal. Classification should therefore specify scientific meaning
and lifetime, not force one GPU kernel per category.

## Comparison

| Concern | Terrarium | ClimaLand | Agrocosm target |
|---|---|---|---|
| Persistent model state | `prognostic` | `Y` | `state` |
| Rate/net change of state | `tendencies` | `dY` | explicit daily increments or conservation updates |
| Scientific daily flux | generally `auxiliary` | cache `p` | `fluxes` |
| Derived limiting factor or geometry | `auxiliary` | cache `p` | `auxiliary` |
| Discrete daily transition marker | generally `auxiliary` | cache `p` | `events` |
| Kernel scratch storage | auxiliary field or time-stepper cache | cache `p` | `workspace` |
| External forcing | `inputs` | drivers/boundary conditions | `inputs` outside crop state |
| Component hierarchy | processes plus `namespaces` | nested component models | crop subcomponents plus a coupling view |
| Output selection | fields selected by user/output writer | diagnostics system reads `Y` or `p` | output registry independent of storage category |

## Recommended Agrocosm design

Agrocosm should keep five distinct lifecycle categories because a daily crop
model benefits from stronger lifecycle and conservation semantics than either
framework's generic auxiliary bucket:

```text
Crop
├── state
│   ├── phenology
│   ├── canopy
│   ├── carbon
│   ├── nitrogen
│   ├── water
│   └── calendar
├── fluxes
│   ├── carbon
│   ├── nitrogen
│   └── water
├── auxiliary
│   ├── canopy
│   ├── photosynthesis
│   └── stress
├── events
│   ├── sowing
│   └── harvest
└── workspace
```

This is a semantic structure, not necessarily five kernel launches or five
physical allocations. Kernels may read and write several categories when that
is demonstrably faster and preserves the declared lifetime rules.

The classification rules are:

1. A field is `state` if tomorrow's solution cannot be reconstructed without
   today's stored value. This includes stocks, phenological memory, event
   history needed by the next day, and management queues that survive a day.
2. A field is `fluxes` if it is an amount or rate transferred during the current
   day and is reset at the daily boundary. Fluxes participate in C/N/water or
   energy balance diagnostics.
3. A field is `auxiliary` if it can be recomputed from current state, current
   inputs, and parameters without using its previous value. An auxiliary may
   be scientifically important and may be selected for output.
4. A field is `events` if it marks a discrete transition that occurred during
   the current day, such as sowing or harvest. Events reset daily but are not
   physical transfers.
5. A field is `workspace` if it exists only to avoid allocation or recomputation
   inside kernels. It is never a restart variable and is not a stable public
   output.
6. Parameters and external inputs are not placed in these lifecycle groups.
7. Output/checkpoint classification is orthogonal: checkpoints save state;
   output may select state, fluxes, or auxiliary values; workspace is excluded.

Because Agrocosm currently uses a one-day discrete transition rather than a
general ODE integrator, `state` still means prognostic state. Its mathematical
form is `x(t + 1 day) = F(x(t), u(t), p)`. We do not need to force daily crop
logic into a continuous tendency API merely to use the same terminology.

## Terrarium coupling path

Coupling is feasible and the two projects have unusually well-aligned goals:
CPU/GPU portability, differentiability, column/grid batching, and modular
process composition. The safest design is a thin adapter rather than making
Agrocosm's core depend directly on Oceananigans or Terrarium.

The proposed mapping is:

- Agrocosm `state` -> Terrarium prognostic fields;
- Agrocosm daily net stock changes -> Terrarium tendencies or a scheduled
  discrete daily update, depending on the Terrarium integrator contract chosen
  for crops;
- Agrocosm `fluxes` and `auxiliary` -> Terrarium auxiliary variables;
- selected Agrocosm workspace -> private/preallocated auxiliary fields only
  when the Terrarium execution path needs them;
- climate and soil boundary data -> Terrarium input variables;
- crop, soil, and management groups -> Terrarium namespaces/components;
- C/N/water/energy exchange -> explicit coupling/boundary flux interfaces.

Two integration modes should remain possible:

1. **Agrocosm-owned column:** Agrocosm advances crop and soil internally while
   Terrarium provides grid, forcing, clock, and orchestration.
2. **Terrarium-owned land column:** Agrocosm supplies crop processes while
   Terrarium owns soil water/energy/carbon. This requires a narrower crop--soil
   exchange contract and is the more valuable long-term architecture.

The second mode should exchange quantities rather than structs. A first version
of the boundary contract should include root-zone water availability, root
uptake by layer, transpiration, canopy interception/throughfall, litter C/N by
layer, mineral-N uptake by layer, surface energy/radiation terms, and management
inputs. This prevents either package from reaching into the other's internal
state layout.

Terrarium is currently early-stage and its interfaces may change. Agrocosm
should therefore stabilize semantic metadata and an adapter protocol before
depending on Terrarium types. A package extension or a small bridge package is
preferable to importing Terrarium in Agrocosm's core. Terrarium is EUPL-1.2,
whereas Agrocosm is MIT; keeping the bridge boundary explicit also makes
licensing and ownership of modifications easier to review.

## Implications for the upcoming audit

The first structural pass now moves the fields into the five lifecycle groups.
A later metadata audit can add a machine-readable registry containing, for
every crop field:

- category and physical units;
- spatial shape (`cell`, `layer × cell`, or `organ × cell`);
- owner process and writer stage;
- daily reset rule;
- restart requirement;
- default output status;
- conservation budget membership;
- expected Terrarium mapping.

The structural hierarchy remains Agrocosm-owned. Terrarium integration is an
optional adapter that maps this public scientific model; it is not a reason to
change Agrocosm's storage semantics or make Terrarium a dependency.
