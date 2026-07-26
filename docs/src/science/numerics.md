# Numerics and conservation

This page documents the ordered transition, backend execution model, scalar
solvers, stability guards, and daily conservation diagnostics.

## Daily transition

For prognostic state ``x_d``, forcing ``u_d``, management ``m_d``, and
parameters ``\vartheta``, one model day is

```math
x_{d+1}=\mathcal{T}(x_d,u_d,m_d;\vartheta),
```

with diagnostics and boundary fluxes

```math
y_d=\mathcal{H}(x_d,u_d,m_d;\vartheta).
```

The transition is ordered rather than a simultaneous algebraic solve. Existing
organic matter decomposes before crop uptake, so newly mineralized nitrogen is
available that day. Harvest residue is routed after existing litter
decomposition and therefore first decomposes the next day. The executable
sequence is listed in [Daily process order](../concepts/daily_processes.md).

## Lifecycle state

Agrocosm separates prognostic variables, daily fluxes, auxiliaries, inputs,
events, workspace, and retained output. A quantity belongs to prognostic state
when its current value affects a later transition. This includes rolling
climate memory, saturation fraction, soil enthalpy, and previous-step thermal
properties; they cannot be reclassified as end-of-day diagnostics merely
because they can be recomputed locally.

## Cell-local backend execution

The reference and accelerator pathways call the same cell-local equations
through KernelAbstractions. One work item owns one horizontal cell and loops
over the fixed vertical soil column where dependencies exist. This avoids
races between layers of the same column and keeps the process implementation
backend-independent.

The batch layout is structure-of-arrays: horizontal cells form the parallel
dimension and soil layers are the small inner dimension. Global runs retain a
fixed compact cell ordering keyed by canonical `cellid`; forcing is streamed
in time blocks while model state stays resident on the selected backend.

## Internal-CO₂ root solve

When water supply constrains conductance, realized internal CO₂ is obtained
from

```math
F(\lambda)=f_g^*(1-\lambda)-A_{dt,mm}(\lambda)=0.
```

Each cell uses bounded bisection over the admissible ``\lambda`` interval. The
iteration count is capped, the bracket never leaves the physical interval, and
the same scalar procedure is used for C3 and C4. This produces deterministic
control flow without allocating a global nonlinear-solver object.

## Exact decay and bounded withdrawals

First-order pool loss is evaluated as

```math
\Delta X=-\operatorname{expm1}(-kR)X,
```

which is more accurate than ``(1-\exp(-kR))X`` for small Float32 rates. Water,
carbon, and nitrogen withdrawals are capped by source storage. Fractions and
environmental response functions are clipped to declared domains, infiltration
uses a bounded number of slugs, and photosynthesis protects its square-root
discriminant:

```math
\Delta_{photo}=\max\left[0,(J_e+J_c)^2-4\theta J_eJ_c\right].
```

Failed crops are terminated before negative living biomass persists.

## Conservation residuals

For conserved quantity ``Q``, the daily residual is

```math
\varepsilon_Q=Q_{d+1}-Q_d-\sum F_{in}+\sum F_{out}.
```

Carbon boundaries include assimilation, harvest export, autotrophic
respiration, and heterotrophic respiration. Nitrogen includes prescribed or
automatic fertilizer, manure, harvest export, leaching, N₂, N₂O, and NH₃.
Water includes precipitation, irrigation, runoff, drainage, interception,
transpiration, and evaporation. Energy includes surface forcing and advective
enthalpy carried by water.

For scale ``S_Q`` the reported relative diagnostic is conceptually

```math
r_Q=\frac{|\varepsilon_Q|}{\max(S_Q,\epsilon)}.
```

Large energy stocks can produce visibly nonzero absolute Float32 residuals
while relative closure remains tight, so absolute and scale-aware residuals
must be inspected together. CPU/GPU equivalence is separate from physical
closure: two backends can agree while sharing a missing boundary term.

## Differentiable boundary

The public one-day `transition_day!` is the intended primal boundary for
future differentiation. Input loading, agricultural warm-up, checkpoint I/O,
streamed output, and reporting remain outside it. Discrete sowing, harvest,
fertilization, clamps, bisection, and crop failure require explicit gradient
policies before the Enzyme CPU path can be considered validated.

## Main safeguards

- non-negative photosynthesis discriminants;
- bounded bisection for ``\lambda``;
- exact exponential decay through `expm1`;
- storage-limited water, carbon, and nitrogen withdrawals;
- explicit bounds on fractions and response functions;
- bounded infiltration iterations;
- crop termination before negative living pools persist.

## Code map

- `src/simulations/simulation_api.jl`
- `src/simulations/runtime_contracts.jl`
- `src/diagnostics/carbon_balance.jl`
- `src/diagnostics/nitrogen_balance.jl`
- `src/diagnostics/water_balance.jl`
- `src/diagnostics/thermal_balance.jl`
- `src/processes/crop/lambda_solver.jl`
- `src/processes/soil/infil_perc.jl`
