# State-classification refactor test plan

This refactor moves restart-critical crop variables into `crop.state`, creates
`crop.auxiliary.root`, and computes top-three-layer root-zone available water
inside the existing transpiration kernel. Run the focused checks first, then
the complete CPU suite. Do not compare results with pre-refactor checkpoints:
the field paths have intentionally changed, while daily numerical behaviour is
expected to remain unchanged.

Run all commands from the repository root.

## Focused CPU checks

```bash
julia --project=. test/processes/initialization/test_initialization.jl
julia --project=. test/processes/crop/test_actual_lai.jl
julia --project=. test/processes/crop/test_crop_lifecycle.jl
julia --project=. test/processes/crop/test_field_lifecycle.jl
julia --project=. test/processes/crop/test_prescribed_phenology.jl
julia --project=. test/processes/crop/test_lambda_water_coupling.jl
julia --project=. test/processes/crop/test_nitrogen_uptake.jl
julia --project=. test/processes/soil/test_soil_water.jl
julia --project=. test/simulations/test_simulation_api.jl
julia --project=. test/simulations/test_cross_year_continuity.jl
julia --project=. test/simulations/test_daily_crop_C3_precision.jl
```

`test_soil_water.jl` includes the new root-zone diagnostic assertion:

```julia
zone_available_water == sum(relative_content[1:3] .* holding_capacity_storage[1:3] .* root_distribution[1:3])
```

## Complete CPU suite

```bash
julia --project=. test/runtests.jl
```

## Focused GPU checks

```bash
julia --project=. test/processes/initialization/test_initialization_gpu.jl
julia --project=. test/processes/crop/test_actual_lai_gpu.jl
julia --project=. test/processes/crop/test_crop_lifecycle_gpu.jl
julia --project=. test/processes/crop/test_daily_process_kernels_gpu.jl
julia --project=. test/processes/crop/test_process_kernels_gpu.jl
julia --project=. test/processes/crop/test_nitrogen_uptake_gpu.jl
julia --project=. test/simulations/test_daily_crop_C3_gpu.jl
julia --project=. test/simulations/test_daily_crop_C3_precision_gpu.jl
```

## Remaining GPU regression scripts

```bash
julia --project=. test/processes/crop/test_fertilizer_gpu.jl
julia --project=. test/processes/crop/test_lambda_solver_c3_gpu.jl
julia --project=. test/processes/crop/test_lambda_solver_c4_gpu.jl
julia --project=. test/processes/crop/test_nitrogen_allocation_gpu.jl
julia --project=. test/processes/crop/test_nitrogen_demand_gpu.jl
julia --project=. test/processes/crop/test_nitrogen_vcmax_limit_gpu.jl
julia --project=. test/processes/crop/test_canopy_snow_cover_gpu.jl
julia --project=. test/processes/soil/test_c_shift_routing_gpu.jl
julia --project=. test/processes/soil/test_litter_routing_gpu.jl
julia --project=. test/processes/soil/test_nitrogen_transform_gpu.jl
julia --project=. test/processes/soil/test_percolation_enthalpy_gpu.jl
julia --project=. test/processes/soil/test_soil_cn_decomposition_gpu.jl
julia --project=. test/processes/soil/test_soil_decomposition_fluxes_gpu.jl
julia --project=. test/processes/soil/test_soil_decomposition_response_gpu.jl
julia --project=. test/processes/soil/test_soil_freeze_thaw_gpu.jl
julia --project=. test/processes/soil/test_soil_process_kernels_gpu.jl
julia --project=. test/processes/soil/test_soil_temperature_gpu.jl
julia --project=. test/processes/soil/test_surface_litter_water_gpu.jl
julia --project=. test/processes/soil/test_untracked_water_enthalpy_gpu.jl
julia --project=. test/processes/soil/test_water_ice_pools_gpu.jl
julia --project=. test/diagnostics/test_nitrogen_balance_gpu.jl
```
