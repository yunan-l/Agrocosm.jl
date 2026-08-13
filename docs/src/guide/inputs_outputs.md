# Inputs and outputs

## Initial-data schema

`initialize_simulation` accepts raw input data with these top-level fields:

- `latitude`
- `crop`: `sdate`, `phu`, `manure`, `fertilizer`, `residuefrac`
- `soilparam`: `soilph`, `w_sat`, `sand`, `clay`, `tdiff_0`, `tdiff_15`, `soildepth`
- `initial_state`: `swc`, `litc`, `fastc`, `slowc`, `litn`, `fastn`, `slown`

Mineral nitrate/ammonium are initialized from slow organic N by default.
Explicit restart pools and post-spin-up routing are available through the
lower-level loader options.

## Climate schema

Raw climate uses daily matrices `(day, cell)`:

- `temp`: air temperature in °C
- `prec`: precipitation in mm day⁻¹
- `swdown`: downward shortwave radiation
- `lwnet`: net longwave radiation
- `windspeed`: optional wind speed; a default is used when absent
- `co2`: annual vector or daily matrix
- `temp_spinup`: temperature history used only to initialize climate memory
- `prec_spinup`: optional precipitation history aligned with `temp_spinup`;
  required to initialize the precipitation climatology before the first
  `:dynamic_sdate` year

## Sowing mode

`initialize_simulation(...; sowing_mode=:prescribed_sdate)` is the default and
uses `crop.sdate` directly. Set `sowing_mode=:dynamic_sdate` to resolve daily
establishment for wheat, maize, rice, or soybean from the climate buffer. This
mode preserves `crop.sdate` as the bootstrap/reference date and never mutates
the source management input. Rain-season timing uses the long-term monthly
precipitation climatology currently available to Agrocosm.

For the global TOML workflows, select the same mode under `[management]`:

```toml
[management]
sowing_mode = "dynamic_sdate"
```

With fixed management, prescribed sowing retains the LPJmL-style fixed
vernalization requirement. Dynamic sowing instead updates it from the annual
climate buffer.

Climate blocks may be passed as `NamedTuple`s, JLD2 paths, or an ordered
vector of blocks. Output time remains continuous across blocks.

`AgrocosmData.prefetch_climate_forcings(reader)` reads the following climate
block on a worker thread while the current block is simulated. Start Julia
with at least two threads for actual overlap.

## Output groups

- `simulation.output.crop`: GPP, NPP, LAI, biomass, yield, Vcmax, plant
  respiration, vegetation C/N, and water deficit. `crop.gpp` is the daily
  gross-assimilation flux in `gC m-2 day-1`.
- `simulation.output.soil`: soil water, soil/litter C and N, and the daily
  ecosystem flux diagnostics below.

### Daily ecosystem flux diagnostics

The daily variables below are model diagnostics; they do not alter any process
calculation.

- `soil.ecosystem_respiration` is modeled **RECO** in `gC m-2 day-1`:
  `crop.respiration + crop.leaf_respiration + soil.heterotrophic_respiration`.
- `soil.heterotrophic_respiration` is the litter and soil microbial component
  of RECO in `gC m-2 day-1`.
- `soil.evapotranspiration` is modeled **ET_total** in `mm day-1`: layered
  crop transpiration, bare-soil evaporation, surface-litter evaporation, and
  canopy-interception evaporation.
- `crop.yield` is emitted annually in `gC m-2 year-1`, at the calendar-year
  boundary. `calendar.harvest_date` and `calendar.harvesting_year` identify
  the harvest event used for that annual record, including a winter crop whose
  sowing and harvest occur in different calendar years.
- `simulation.output.climate`: selected processed climate variables.
- `simulation.output.calendar`: sowing/harvest events and harvest dates.

Soil and climate output coverage is still being expanded. For scientific
state inspection, use the lifecycle tree rather than assuming every internal
field has a time-series output.

For global runs, attach an `OutputStream`. Equal-sized output blocks reuse the
same backend buffers; emitted rows are reduced or written before the next
block and are not retained in model memory. NetCDF stream variables include
their units and long names from the output-variable schema.

The numerical state inventory is available through `state_schema`, whose
entries identify the lifecycle role and dimensions of each array.
