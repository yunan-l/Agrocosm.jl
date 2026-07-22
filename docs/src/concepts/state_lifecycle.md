# State lifecycle

Agrocosm separates scientific process organization from numerical variable
lifetime. Both crop and soil variables use the same canonical lifecycle tree:

```text
ModelState
├── prognostic
├── fluxes
├── auxiliary
├── inputs
├── events
├── workspace
└── output
```

## Categories

| Group | Contract |
|---|---|
| `prognostic` | Cross-day memory required to determine the next transition. |
| `fluxes` | Transfers produced for the current day and overwritten by owner processes. |
| `auxiliary` | Algebraic or diagnostic fields derived during the current process chain. |
| `inputs` | External forcing, static geometry, management input, or retained configuration. |
| `events` | One-day discrete sowing and harvest indicators. |
| `workspace` | Preallocated numerical scratch with no scientific meaning. |
| `output` | User-facing copies recorded after a completed daily transition. |

Examples:

```julia
state.prognostic.crop.carbon.leaf
state.prognostic.soil.water.storage
state.fluxes.crop.carbon.npp
state.fluxes.soil.water.percolation
state.auxiliary.crop.photosynthesis.vcmax
state.inputs.weather.temp
state.events.crop.harvest
```

An auxiliary field must not silently become cross-day memory. Under the
current ordering, snow height/fraction and prior thermal coefficients affect
the following transition and are therefore classified as prognostic state.

Only prognostic variables and restart-relevant inputs belong to the scientific
checkpoint. Fluxes, recomputable auxiliaries, and workspaces do not.

The exhaustive field inventory, including units, is maintained in the source
repository's `docs/variable_inventory.md`.
