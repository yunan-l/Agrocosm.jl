# Validation and limitations

## What is tested

- CPU process and integration tests for crop, soil water, C/N, heat, outputs,
  checkpoints, and daily ordering.
- Dedicated CUDA process and C3 end-to-end equivalence scripts.
- `Float32` and `Float64` paths.
- Daily water, carbon, nitrogen, thermal-energy, and percolation-energy ledgers.
- Interrupted/checkpointed trajectories against uninterrupted simulations.

The full local CPU `Pkg.test()` suite passes. The three-day C3/C4 lifecycle
migration was compared against the previous runtime entry point across all
arrays before that entry was removed (`1032/1032` exactly equal).

## Scientific interpretation

Agrocosm is research software. Passing conservation and implementation tests
does not establish universal agronomic validity. Parameter sets must be
evaluated for the crop, cultivar, management system, soil, climate, spatial
scale, and question of interest.

Current limitations include:

- no equilibrium soil/ecosystem spin-up workflow; the available finite
  agricultural warm-up does not equilibrate the slow SOC pool;
- the multi-CFT patch runner is implemented, but its complete 24-patch global
  production matrix still needs server evidence;
- no production Penman–Monteith/Medlyn alternative;
- simplified frozen-soil infiltration and heat transport;
- incomplete soil/climate time-series output coverage;
- no public differentiable one-day transition yet;
- no broad multi-site or global validation protocol.

Initial stock drift in long simulations without spin-up should not be
interpreted as equilibrium behaviour.
