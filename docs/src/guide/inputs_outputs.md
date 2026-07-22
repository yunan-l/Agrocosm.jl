# Inputs and outputs

## Initial-data schema

`initialize_simulation` accepts raw input data with these top-level fields:

- `latitude`
- `crop`: `sdate`, `phu`, `manure`, `fertilizer`, `residuefrac`
- `soilparam`: `soilph`, `w_sat`, `sand`, `clay`, `tdiff_0`, `tdiff_15`, `soildepth`
- `initialLPJmL.u0`: `swc`, `litc`, `fastc`, `slowc`, `litn`, `fastn`, `slown`

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

Climate blocks may be passed as `NamedTuple`s, JLD2 paths, or an ordered
vector of blocks. Output time remains continuous across blocks.

## Output groups

- `simulation.output.crop`: GPP, NPP, LAI, biomass, yield, Vcmax, respiration,
  vegetation C/N, and water deficit.
- `simulation.output.soil`: soil water, soil/litter C and N, ecosystem and
  heterotrophic respiration, and evapotranspiration.
- `simulation.output.climate`: selected processed climate variables.
- `simulation.output.calendar`: sowing/harvest events and harvest dates.

Soil and climate output coverage is still being expanded. For scientific
state inspection, use the lifecycle tree rather than assuming every internal
field has a time-series output.
