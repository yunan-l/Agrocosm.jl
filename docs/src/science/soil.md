# Soil processes

This page collects the soil water, thermal, carbon, and nitrogen equations.

## Water

The soil column is currently five layers deep. Rain and snowmelt enter the top
boundary after interception. Infiltration is processed in bounded slugs for
numerical stability. Layer water is capped by storage limits, and excess water
becomes percolation or runoff. Drainage exports the source-layer water flux.

Full irrigation, when enabled, restores the rooted soil to field capacity every
day without an external water-availability constraint.

## Heat

Soil temperature is not evolved independently from the other soil state.
Instead, the code reconstructs temperature and phase partition from enthalpy
and layer thermal properties. Liquid water carries enthalpy between layers, and
the model recomputes heat capacity and conductivity after each water transfer.

The soil thermal state is therefore a lifecycle variable, not a pure daily
diagnostic.

## Decomposition response

The soil decomposition response combines a temperature factor and a moisture
polynomial. Surface litter uses surface-litter temperature and wetness, while
incorporated and below-ground litter use the topsoil response.

## Carbon pools

Litter, fast, and slow pools use exact daily first-order decay,

```math
\Delta C=(1-e^{-kR})C,
```

implemented with `expm1` for numerical stability. Litter loss is split between
atmospheric respiration and retained material, and heterotrophic respiration is
the sum of the atmospheric share and the decomposed fast/slow pools.

Tillage modifies topsoil hydraulic properties and routes litter differently from
the bioturbation pathway.

## Nitrogen pools

Organic litter, fast, and slow nitrogen follow the same response and routing as
their carbon counterparts. Gross mineralization enters ammonium. If litter carbon
requires more nitrogen than decomposition supplies, the model immobilizes
mineral nitrogen back into organic pools.

Nitrification, denitrification, and NH₃ volatilization are then applied in the
daily soil-N sequence. Volatilization acts on the whole top-layer ammonium
pool, which is why the cumulative loss can remain large when mineralization and
fertilization keep replenishing ammonium.

## Code map

- `src/processes/soil/soil_water.jl`
- `src/processes/soil/infil_perc.jl`
- `src/processes/soil/evaporation.jl`
- `src/processes/soil/soil_temp.jl`
- `src/processes/soil/soil_response.jl`
- `src/processes/soil/soil_carbon.jl`
- `src/processes/soil/nitrogen_transform.jl`
- `src/processes/soil/litter_routing.jl`
- `src/processes/soil/tillage.jl`
