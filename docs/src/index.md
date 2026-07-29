# Agrocosm.jl

**Agrocosm.jl is a process-based, GPU-capable Julia model of crop–soil carbon,
water, nitrogen, and energy dynamics.** It is designed for transparent daily
simulation from individual sites to batches of independent grid cells.

The current model provides C3 and C4 crop pathways, phenology and management,
five-layer soil hydrology and heat transport, snow and freeze–thaw processes,
and coupled soil carbon–nitrogen transformations. Carbon, nitrogen, water, and
thermal balance ledgers make the implementation auditable.

Agrocosm uses LPJmL crop and soil processes as an important scientific
reference, but is an independent Julia implementation rather than a
line-by-line translation.

## Why Agrocosm?

- One process implementation runs on CPU and GPU backends.
- Both `Float32` and `Float64` simulations are supported.
- Process configuration is separated from the numerical State variables.
- File checkpoints can resume a simulation at completed daily boundaries.
- Conservation diagnostics are built into end-to-end simulations.

## Current maturity

The model can run up to 12 LPJmL-compatible CFTs under rainfed and full
irrigation modes (up to 24 independent crop-soil patches). The multi-CFT CPU
runner writes one unweighted yield band per selected CFT and water system.
The model remains under active research development: crop rotations,
Penman–Monteith/Medlyn alternatives, and end-to-end automatic differentiation
are not yet production features.

Start with [Getting started](@ref), then read [Model overview](@ref) and
[State variables](@ref) before extending a process.

```@contents
Pages = [
    "getting_started.md",
    "concepts/overview.md",
    "concepts/state_lifecycle.md",
    "concepts/daily_processes.md",
]
Depth = 2
```
