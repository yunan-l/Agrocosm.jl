# API reference

## Simulation interface

```@docs
initialize_simulation
run_simulation!
simulation_summary
save_checkpoint
restore_checkpoint!
CropSimulation
ProcessModules
ModelState
```

## Data preparation

```@docs
InitialDataLoader
ClimateDataLoader
```

## Main process entry points

```@docs
daily_crop_C3!
daily_crop_C4!
pedotransfer!
soil_temperature!
soil_cn_decomposition!
photosynthesis_C3!
photosynthesis_C4!
transpiration!
crop_carbon!
crop_nitrogen!
```
