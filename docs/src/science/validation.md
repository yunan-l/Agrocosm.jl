# Validation and limitations

## What is tested

- CPU process and integration tests for crop, soil water, C/N, heat, outputs,
  checkpoints, and daily ordering.
- Dedicated CUDA process and C3 end-to-end equivalence scripts.
- `Float32` and `Float64` paths.
- Daily water, carbon, nitrogen, thermal-energy, and percolation-energy ledgers.
- Interrupted/checkpointed trajectories against uninterrupted simulations.

The current post-refactor CPU suite contains 1686 passing assertions. The
three-day C3/C4 lifecycle migration was compared against the previous runtime
entry point across all arrays before that entry was removed (`1032/1032`
exactly equal).

## Scientific interpretation

Agrocosm is research software. Passing conservation and implementation tests
does not establish universal agronomic validity. Parameter sets must be
evaluated for the crop, cultivar, management system, soil, climate, spatial
scale, and question of interest.

Current limitations include:

- no completed soil/ecosystem spin-up workflow;
- no multi-crop competition or rotation framework;
- no production Penman–Monteith/Medlyn alternative;
- simplified frozen-soil infiltration and heat transport;
- incomplete soil/climate time-series output coverage;
- no public differentiable one-day transition yet;
- no broad multi-site or global validation protocol.

Initial stock drift in long simulations without spin-up should not be
interpreted as equilibrium behaviour.
