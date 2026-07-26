# Model equations processes

This section is the entry point for Agrocosm's science documentation. The
details now live in smaller pages grouped by process domain:

- [Crop processes](crop.md)
- [Soil processes](soil.md)
- [Climate and surface processes](climate_surface.md)
- [Numerics and conservation](numerics.md)
- [Initialization and warm-up](initialization_warmup.md)

Agrocosm advances a batch of independent crop–soil columns at a daily time
step. The implementation is LPJmL-informed, but the equations below describe
Agrocosm's current code path rather than an abstract textbook model.

## Scope

The current production configuration assumes:

- one crop PFT per stand and grid cell;
- a five-layer soil column;
- daily forcing on a strict 365-day, no-leap calendar;
- rainfed conditions or unconstrained full irrigation;
- C3 or C4 photosynthesis;
- coupled crop carbon, crop nitrogen, soil water, soil heat, and soil C/N;
- no natural-vegetation competition or lateral exchange between grid cells;
- no dedicated frozen-soil infiltration impedance or frozen-soil heat advection.

## State and outputs

Agrocosm separates prognostic state, fluxes, auxiliaries, inputs, and discrete
events. That split matters for restart behavior, backend portability, and the
future differentiable transition API.

In one daily step, the model can be written abstractly as

```math
x_{d+1}=\mathcal{T}(x_d,u_d,m_d;\vartheta),\qquad
y_d=\mathcal{H}(x_d,u_d,m_d;\vartheta),
```

where `x` is the evolving state, `u` the daily forcing, `m` the management
input, and `\vartheta` the model parameters.

## Daily order

The executable order is documented in [Daily process order](../concepts/daily_processes.md).
At a high level, each day updates:

1. climate memory and management;
2. snow, albedo, radiation, and PET;
3. soil water, soil heat, and soil C/N;
4. crop phenology and harvest logic;
5. canopy interception, photosynthesis, transpiration, and nitrogen uptake;
6. respiration, allocation, and biomass pools;
7. evapotranspiration, fertilizer/manure, and post-crop nitrogen losses.

Order matters. Existing organic matter can mineralize and become available for
same-day crop uptake. New harvest residues are routed after that day's current
litter decomposition, so they only decompose on the following day.

## Units

Unless stated otherwise, carbon and nitrogen fluxes use `gC m⁻² day⁻¹` and
`gN m⁻² day⁻¹`, water uses `mm day⁻¹`, and state variables are expressed per
unit ground area.

## Where to go next

If you want the actual equations, start with [Crop processes](crop.md) and
[Soil processes](soil.md). If you want the backend and closure assumptions,
read [Numerics and conservation](numerics.md). If you want the global-input and
warm-up assumptions, read [Initialization and warm-up](initialization_warmup.md).
