# Soil C–N decomposition and daily-flow audit

This note records the LPJmL behavior used by Agrocosm's single-crop,
non-methane soil C–N pathway. The local LPJmL reference is commit
`572e2b906ac2c55b2ee6661a93e4633b126254e4`.

## LPJmL source basis

- `src/lpj/update_daily_cell.c`
  - sowing and stand replacement occur before the stand loop;
  - daily bioturbation is applied before soil thermal and decomposition work;
  - `daily_littersom()` is called before `daily_stand()`;
  - `denitrification()` and `volatilization()` are called after `daily_stand()`.
- `src/soil/daily_littersom.c`
  - selects the non-methane `littersom_nomethane()` path when methane is off.
- `src/soil/littersom_nomethane.c`
  - applies the same pool-specific fractional decay to C and N;
  - mineralizes decomposed SOM N into NH4;
  - partitions non-respired litter C and N into fast/slow pools with `c_shift`;
  - mineralizes the `atmfrac` share of decomposed litter N;
  - immobilizes mineral N when decomposed litter C:N exceeds the soil target;
  - nitrifies after mineralization/immobilization;
  - does not activate the commented TODO for N-limited SOM decay.

## Agrocosm daily order after this audit

```text
1. cultivation, prescribed inputs, tillage, bioturbation
2. surface-litter properties and five-layer soil temperature
3. pre-crop soil C–N stage
   a. shared temperature/moisture decomposition response
   b. existing litter and fast/slow SOM C and N decomposition
   c. litter-to-SOM transfer through fast/slow c_shift
   d. SOM/litter N mineralization and litter-N immobilization
   e. nitrification
4. phenology and harvest
5. route new harvest C and N residues; do not decompose them today
6. infiltration, crop water/C/N processes, evaporation and soil-water removal
7. post-crop denitrification and NH3 volatilization
8. C, N, water and energy balance diagnostics
```

The split is implemented by:

- `soil_cn_decomposition!`: pre-crop decomposition, mineralization,
  immobilization and nitrification;
- `route_harvest_residues!`: same-day conservative C/N residue routing;
- `post_crop_nitrogen_losses!`: post-crop denitrification and volatilization.

The older `soil_carbon!`, `soil_nitrogen!`, and `nitrogen_transform!` entry
points remain available for isolated-process compatibility tests.

## Conservation and coupling invariants

For each litter or SOM pool with positive C and N, C and N use the same decay
fraction. During the pre-crop stage:

```text
initial C = final C + heterotrophic respiration
initial N = final organic N + final mineral N + nitrification N2O
```

Mineralization, immobilization, nitrification, denitrification and
volatilization are reported separately. Immobilization may not make NH4 or
NO3 negative and is limited by the mineral-N supply in each layer.

## Deliberately retained differences

- Agrocosm currently represents one crop/litter class per simulated cell,
  rather than LPJmL's list of PFT-specific litter items. It therefore uses one
  fast and one slow `c_shift` profile for the active crop case.
- Leaf-like surface, incorporated and root litter are represented, but woody
  fuel-size classes are not needed for the current annual crop scope.
- The current pathway is the LPJmL non-methane pathway. O2/CH4 transport,
  adaptive sub-daily methane timesteps, wetlands and rice methane remain out
  of scope until the rice/wetland stage.
- Atmospheric N deposition and biological N fixation are not yet part of this
  single wheat case.
- Agrocosm uses `-expm1(-x)` for exponential decay to reduce Float32
  cancellation; this is mathematically equivalent to LPJmL's `1-exp(-x)`.

## Regression tests

- `test/processes/soil/test_soil_cn_decomposition.jl`
  - identical C/N decay fractions;
  - carbon and nitrogen conservation;
  - N-rich net mineralization behavior;
  - N-poor immobilization and mineral-supply limitation;
  - separation of pre-crop and post-crop nitrogen stages.
- `test/processes/soil/test_soil_cn_decomposition_gpu.jl`
  - CPU/CUDA agreement for all C/N pools and transformation fluxes.
- `notebooks/wheat_rainfed_api.ipynb`
  - 20-year high-level end-to-end regression.
