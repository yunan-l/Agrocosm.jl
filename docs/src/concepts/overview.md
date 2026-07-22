# Model overview

Agrocosm advances a batch of independent crop–soil columns at a daily time
step. Array leaves use a structure-of-arrays layout: cell fields have shape
`(cell)`, soil fields `(layer, cell)`, and output fields `(day, cell)`.

## Crop processes

- C3 and C4 photosynthesis with water-limited internal CO₂ ratio.
- Canopy radiation interception, phenology, LAI development, and senescence.
- Maintenance and growth respiration and organ carbon allocation.
- Plant nitrogen demand, separate nitrate/ammonium uptake, and organ allocation.
- Cultivation, split fertilizer/manure application, harvest, and residue return.

## Soil processes

- Five-layer water storage, infiltration, runoff, percolation, evaporation, and
  root water removal.
- Snow accumulation/melt and surface litter interception.
- Soil temperature, enthalpy, freeze–thaw partitioning, and water-carried heat.
- Litter, fast, and slow carbon/nitrogen pools.
- Decomposition, mineralization, immobilization, nitrification,
  denitrification, ammonia volatilization, and mineral-N leaching.
- Tillage effects on litter routing and topsoil hydraulic properties.

## Process configuration and numerical state

`ProcessModules` owns process choices and parameters. `ModelState` owns all
evolving and diagnostic arrays. Process wrappers select lifecycle-scoped
arrays from `ModelState` and pass explicit leaves to backend kernels. This
separation is the foundation for interchangeable processes and future
differentiable transitions.

The current production interface advances a range of days. A dedicated
one-day differentiable transition and Enzyme integration are planned but not
yet part of the public API.
