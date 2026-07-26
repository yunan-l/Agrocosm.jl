# Numerics and conservation

This page explains the numerical choices that sit underneath the process
equations.

## Cell-local execution

Agrocosm launches cell-local kernels through KernelAbstractions. One work item
owns one horizontal grid cell and loops over the vertical soil column where
needed. This keeps CPU and GPU behavior aligned while avoiding cross-cell races.

## Daily coupling

The full daily transition is ordered and stateful. Some quantities must be
computed before others for the model to match LPJmL-style behavior. The main
example is nitrogen: mineralization happens before crop uptake, and harvest
residue is routed after decomposition of the current-day pools.

## Scalar solvers

The internal-CO₂ ratio `\lambda` is solved with bounded bisection. The loop has
a fixed maximum iteration count and is run independently in each cell.

## Stability guards

The implementation relies on a small set of recurring safeguards:

- non-negative discriminants in photosynthesis;
- bounded root finding for `\lambda`;
- exact decay via `expm1`;
- explicit clipping of fractions and environmental response functions;
- storage-limited extraction for water, carbon, and nitrogen;
- bounded infiltration slugs for soil water.

## Conservation diagnostics

Daily carbon, nitrogen, water, and energy balances are tracked separately from
the main state updates. A residual is useful for diagnosing a missing boundary
term or a backend mismatch, but it is not by itself proof of a scientific
problem. Absolute and scale-aware residuals should be read together.

## What is not here yet

- no public one-day differentiable transition API yet;
- no Enzyme integration in the production path;
- no global multi-crop partitioning layer;
- no frozen-soil advective heat transport or impedance model.
