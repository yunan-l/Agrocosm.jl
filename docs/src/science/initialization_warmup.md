# Initialization and warm-up

This page describes the assumptions behind the current input-data path and the
finite agricultural warm-up.

## Input data

Global or regional runs depend on external climate, management, and soil input
files. The current data-preparation workflow extracts the relevant crop PFT,
management bands, climate years, and CO₂ series from larger source datasets so
the model can run on a practical subset without loading the full global
archives into memory.

## Soil initialization

The current production setup initializes soil carbon and nitrogen from HWSD-
based inputs rather than from a full ecosystem spin-up. That keeps Agrocosm
decoupled from LPJmL for the input stage and avoids re-running the source model
whenever the spatial domain changes.

## Agricultural warm-up

Because the model does not yet have a full equilibrium spin-up workflow, it
uses a finite agricultural warm-up before formal production runs. The purpose
is to give litter and fast pools some history and reduce the abruptness of the
first simulated year. It is not a substitute for a full soil-C equilibrium
spin-up.

## Code map

- `lib/AgrocosmData/scripts/prepare_global_wheat_subset.jl`
- `docs/src/guide/global_wheat_subset.md`
- `docs/src/guide/hwsd_initialization.md`
- `src/simulations/agricultural_warmup.jl`
