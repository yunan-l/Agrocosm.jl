# Process contracts

This page records the runtime contract of Agrocosm's core daily processes.
It complements the exhaustive variable inventory: the inventory defines every
field, while this page identifies process owner, inputs, outputs, units, and
event phase.

## Contract conventions

- **Reads** may include preceding state, forcing, static inputs, and parameters.
- **Writes** are owned by the named process for the current day. Flux and
  auxiliary outputs must be overwritten on every relevant path.
- **State** is cross-day and checkpointed. **Flux** is a daily transfer.
  **Auxiliary** is recomputable within the daily transition.
- Carbon is `gC m⁻²`, nitrogen is `gN m⁻²`, water is `mm`, and energy is
  `J m⁻²`, unless stated otherwise.

No process may write a field owned by another process in the same phase. This
keeps each cell independent inside a KernelAbstractions launch and makes the
daily order auditable on CPU and GPU.

## Daily phases

| Phase | Owner processes | Discrete event? | Primary result |
|---|---|---:|---|
| Climate and management | climate history, cultivation, fertilizer/manure | cultivation | Current forcing and stand-management state |
| Surface and thermal preparation | tillage, albedo, snow, pedotransfer, litter properties, soil heat | no | Hydraulic/thermal state for this day |
| Soil C/N turnover | decomposition | no | Organic-pool transfers and mineral N availability |
| Crop calendar | phenology, harvest, residue routing | harvest | Crop phase state and future litter input |
| Water supply | interception, infiltration, percolation | no | Soil water after supply and drainage |
| Canopy assimilation | APAR, temperature stress, photosynthesis, transpiration, `λ` solve | no | GPP and water-limited conductance |
| Crop and soil losses | carbon/N allocation, evaporation, transpiration removal, N losses | crop failure | Updated crop/soil stocks and daily fluxes |
| Diagnostics and output | balances, output writer | no | Completed daily record |

Sowing/cultivation, harvest, and failed-crop termination are explicit one-day
events. They modify state only in their own phase.

## Climate, surface, and soil physics

| Process family | Reads | Owned writes | Lifecycle role | Units | Backend |
|---|---|---|---|---|---|
| Climate history | Daily temperature and preceding `ClimBuf` state | Rolling temperature statistics and vernalization requirement | prognostic climate | °C, day-equivalent | CPU/GPU |
| Snow and surface energy | Precipitation, air temperature, prior snow state | Snow pack, melt/sublimation/runoff fluxes, snow height/fraction | prognostic soil snow; fluxes | mm, m, mm day⁻¹ | CPU/GPU |
| Pedotransfer and litter properties | Texture, tillage state, water state, litter state | Hydraulic auxiliaries; litter cover, depth, and physical properties | auxiliary soil water; prognostic surface litter | mm, 0–1, m | CPU/GPU |
| Soil heat | Air/climate-history temperature, snow/litter properties, prior thermal state | Layer temperature, enthalpy, frozen fraction, thermal fluxes | prognostic soil thermal; fluxes | °C, J m⁻³, J m⁻² day⁻¹ | CPU/GPU |
| Infiltration and percolation | Precipitation, snowmelt, interception, hydraulic auxiliaries, water state | Water storage and infiltration/runoff/percolation fluxes | prognostic soil water; fluxes | mm, mm day⁻¹ | CPU/GPU |
| Evaporation and transpiration removal | Potential evaporation, canopy/litter and soil-water state, transpiration demand | Water storage, evaporation, transpiration-layer fluxes | prognostic soil water; fluxes | mm, mm day⁻¹ | CPU/GPU |

`saturation_fraction` and thermal material coefficients are prognostic because
the next hydraulic or thermal solve reads their previous values. Other
hydraulic and surface-response fields are daily auxiliaries.

## Soil carbon and nitrogen

| Process family | Reads | Owned writes | Lifecycle role | Units | Backend |
|---|---|---|---|---|---|
| Decomposition | Litter, fast/slow C and N pools; temperature/moisture responses; fixed routing shifts | Organic pools, mineral N, decomposition, mineralization, immobilization, heterotrophic respiration | prognostic soil C/N; fluxes and auxiliaries | gC or gN m⁻²; gC or gN m⁻² day⁻¹ | CPU/GPU |
| N transformations and losses | Mineral N, water/thermal state, air temperature, wind | Nitrate/ammonium pools; nitrification, denitrification, volatilization, leaching | prognostic soil N; fluxes | gN m⁻²; gN m⁻² day⁻¹ | CPU/GPU |
| Residue routing | Crop harvest/failure residue fluxes and tillage routing | Litter C/N pools and routing fluxes | prognostic soil C/N; fluxes | gC or gN m⁻²; gC or gN m⁻² day⁻¹ | CPU/GPU |

Decomposition precedes crop uptake, so mineralized N is available on the same
day. Harvest and failed-crop residues enter after decomposition and become
eligible on the following day.

## Crop calendar, canopy, carbon, and nitrogen

| Process family | Reads | Owned writes | Lifecycle role | Units | Backend |
|---|---|---|---|---|---|
| Cultivation and fertilization | Sowing date, PHU/winter type, management, crop/soil state | Sowing event, establishment state, pending inputs, mineral-N and manure-litter additions | events and prognostic crop/soil; fluxes | day of year, gN m⁻² | CPU/GPU |
| Phenology and harvest | Climate history, temperature, day length, phenology state | Heat/vernalization sums, growing/senescence/harvest state, harvest event and exports | prognostic crop; events and fluxes | °C day, day equivalent, gC or gN m⁻² day⁻¹ | CPU/GPU |
| Canopy radiation | LAI, crop/soil/snow albedo, shortwave forcing, day length | Albedo, FPAR, APAR, wet-canopy and interception fluxes | auxiliary crop; fluxes | 0–1, J m⁻² day⁻¹, mm day⁻¹ | CPU/GPU |
| Photosynthesis and conductance | APAR, day length, temperature, CO₂, CFT traits, crop state | `Vcmax`, temperature stress, `λ`, gross/net assimilation, water-limited assimilation | auxiliary crop; fluxes | gC m⁻² day⁻¹, mm day⁻¹, 0–1 | CPU/GPU |
| Crop carbon | Assimilation, respiration drivers, organ C state, soil temperature | Organ/total crop C, NPP, respiration, yield/export fluxes | prognostic crop; fluxes | gC m⁻²; gC m⁻² day⁻¹ | CPU/GPU |
| Crop nitrogen | Soil mineral N, CFT C:N traits, realized `Vcmax`, crop C/N state | Organ/total N, demand/deficit, uptake, automatic-fertilizer flux | prognostic crop; fluxes | gN m⁻²; gN m⁻² day⁻¹ | CPU/GPU |
| Crop failure | Negative-biomass condition, crop/soil state, residue fraction | Failed-crop event, residue fluxes, cleared crop state | events and prognostic crop; fluxes | 0/1, gC or gN m⁻² day⁻¹ | CPU/GPU |

The C3/C4 pathway is selected before kernel launch. Both pathways satisfy the
same state and flux contract; only their biochemical kernels differ.

## Runtime composition and data boundary

Agrocosm uses its own small composition layer rather than an external process
framework:

| Interface | Responsibility | Must not own |
|---|---|---|
| `ProcessModules` | Immutable CFT traits, global parameters, and the C3/C4 pathway selected before daily execution | Numerical arrays or cross-day state |
| `ModelState` | Prognostic, flux, auxiliary, input, event, workspace, and output arrays grouped by lifecycle | Parameter selection or file I/O |
| `ExecutionContext` | Precision, CPU/accelerator array transfer, and stable active-cell identifiers | Scientific process choices |
| `SimulationConfiguration` | Immutable assembly choices: precision, backend, active domain, run duration, and management/process switches | Parameters, arrays, or dataset decoding |
| `CropSimulation` | Assembles the four contracts and advances the audited daily transition | Dataset decoding or preprocessing |
| `AgrocosmData` | NetCDF/database access and backend-neutral climate, management, and soil-input contracts | Warm-up, state transition, or process algorithms |

The resulting boundary is deliberate: file I/O and bounded forcing blocks are
prepared on the CPU, then the selected arrays are transferred once at the
execution boundary. All numerical process work after that point follows the
same CPU/GPU implementation. A future optional adapter may map these explicit
contracts to another model, but no external framework is a dependency of the
Agrocosm runtime.

## Verification contract

At the end of every day, optional water, carbon, nitrogen, thermal, and
percolation-energy ledgers read completed state and owner fluxes. They never
alter process state. Output recording likewise observes, but does not feed
back into, the daily transition.

Changes to a process must update this page when they alter a field owner,
cross-phase input, unit, event timing, or checkpoint requirement. Numerical
equation changes belong in the corresponding model-process documentation.

## Field-level backend equivalence

The CUDA end-to-end test recursively compares every numerical leaf below
`ModelState` after a 730-day C3 run. Float leaves use pointwise tolerances and
integer/boolean leaves, including calendar events, must match exactly. This is
stronger than comparing final yield alone: it covers prognostic state, daily
fluxes, auxiliary fields, inputs, event flags, workspaces, and output arrays.
