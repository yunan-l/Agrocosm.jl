# Climate and surface processes

This page covers the surface and climate-side pieces that feed the crop and
soil kernels.

## Solar geometry and PET

Day length follows the standard astronomical relationship used by LPJmL-style
models. Daily PAR and equilibrium evapotranspiration are diagnosed from short-
wave and long-wave radiation, temperature, and day length.

## Snow

If daily temperature is below the snowfall threshold, precipitation is routed
to the snowpack. Snow is capped at a maximum pack depth, a small fixed
sublimation loss is applied when pack is present, and positive temperature can
melt the pack through the snow-skin energy balance. Melt water is added to the
day's liquid precipitation before interception and infiltration.

## Surface albedo

Effective albedo combines crop canopy, surface litter, bare soil, and snow.
Snow replaces the underlying leaf/litter albedo when the snowpack is present.
When the crop is inactive, the model uses the bare soil/snow mixture only.

## Climate memory

The climate buffer stores short rolling histories used by phenology and related
processes. These buffers are updated on the daily loop and are part of the
numerical state lifecycle.

## Code map

- `src/processes/crop/radiation.jl`
- `src/processes/climate/snow.jl`
- `src/processes/climate/climbuf.jl`
- `src/processes/climate/temp_stress.jl`
- `src/processes/climate/readclimate.jl`
