# API reference

```@meta
CurrentModule = Agrocosm
```

## Managed-crop model

The top-level model and the crop vegetation that slots into a Terrarium `LandModel`.

```@autodocs
Modules = [Agrocosm]
Pages = ["crop/crop_model.jl", "crop/vegetation.jl"]
```

## Crop management

Discrete sowing/harvest events (as Oceananigans callbacks) and continuous fertilizer application.

```@autodocs
Modules = [Agrocosm]
Pages = ["crop/management.jl"]
```

## Crop physiology

C3/C4 photosynthesis, stomatal conductance, phenology, carbon and nitrogen pools, respiration,
allocation, root distribution, and plant-available water.

```@autodocs
Modules = [Agrocosm]
Pages = [
    "crop/photosynthesis.jl",
    "crop/stomatal_conductance.jl",
    "crop/phenology.jl",
    "crop/phenology_dynamics.jl",
    "crop/carbon.jl",
    "crop/carbon_allocation.jl",
    "crop/carbon_dynamics.jl",
    "crop/growth_respiration.jl",
    "crop/maintenance_respiration.jl",
    "crop/nitrogen.jl",
    "crop/nitrogen_allocation.jl",
    "crop/nitrogen_demand.jl",
    "crop/nitrogen_uptake.jl",
    "crop/nitrogen_limitation.jl",
    "crop/harvest_index.jl",
    "crop/root_distribution.jl",
    "crop/plant_available_water.jl",
]
```

## Soil carbon–nitrogen biogeochemistry

The crop soil biogeochemistry and its component processes.

```@autodocs
Modules = [Agrocosm]
Pages = [
    "crop/soil_biogeochemistry.jl",
    "crop/soil_carbon.jl",
    "crop/soil_decomposition_response.jl",
    "crop/nitrification.jl",
    "crop/denitrification.jl",
    "crop/mineralization.jl",
    "crop/volatilization.jl",
]
```

## Crop parameters and the CFT registry

The 12 LPJmL crop functional types and the model parameter sets.

```@autodocs
Modules = [Agrocosm]
Pages = ["parameters/pft.jl", "parameters/default_params.jl", "crop/cft_presets.jl"]
```

## Numerics

```@autodocs
Modules = [Agrocosm]
Pages = ["numerics/lpj_bisect.jl"]
```
