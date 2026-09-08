# Validation and limitations

## What is tested

- CPU process and integration tests for crop, soil water, C/N, heat, outputs,
  checkpoints, and daily ordering.
- Dedicated CUDA process and C3 end-to-end equivalence scripts.
- `Float32` and `Float64` paths.
- Daily water, carbon, nitrogen, thermal-energy, and percolation-energy ledgers.
- Interrupted/checkpointed trajectories against uninterrupted simulations.

Regression tests cover state-lifecycle consistency and numerical equivalence.
Test results apply to the tested revision, inputs, precision and backend;
experiment-specific acceptance records are maintained separately.

## Scientific interpretation

Agrocosm is research software. Passing conservation and implementation tests
does not establish universal agronomic validity. Parameter sets must be
evaluated for the crop, cultivar, management system, soil, climate, spatial
scale, and question of interest.

Current limitations include:

- no equilibrium soil/ecosystem spin-up workflow; the available finite
  agricultural warm-up does not equilibrate the slow SOC pool;
- backend and scientific validation must cover the selected CFT, water-system
  and input configuration;
- no production Penman–Monteith/Medlyn alternative;
- simplified frozen-soil infiltration and heat transport;
- incomplete soil/climate time-series output coverage;
- fixed-event weather derivatives are local model sensitivities, not smooth
  derivatives across discrete sowing, harvest or crop-failure changes;
- full-model weather AD is CPU-only; GPU forcing compatibility is not full
  GPU differentiation support;
- no broad multi-site or global validation protocol.

Initial stock drift in long simulations without spin-up should not be
interpreted as equilibrium behaviour.
