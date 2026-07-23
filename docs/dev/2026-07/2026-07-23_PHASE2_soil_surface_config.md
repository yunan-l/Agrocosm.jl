# Phase 2 — Terrarium soil & surface configuration for the crop model

> Companion to `2026-07-23_PLAN_terrarium_migration.md`. This note records how Agrocosm's legacy
> physical soil/surface/climate code maps onto Terrarium's configurable processes, so the crop
> `LandModel` (Phase 5) can be assembled from Terrarium components and the crop-physics port
> (Phase 3) knows the exact coupling interface. No legacy files are deleted here (retention policy).

Base revision when written: see `git log` for the Phase 2 commit.

## How Terrarium assembles a land column

`LandModel(grid; soil, vegetation, ...)` (`Terrarium/src/models/coupled/land_model.jl:10`) wires
seven `@component` sub-models: `vegetation`, `soil`, `surface_energy_balance`, `surface_hydrology`,
`atmosphere`, `constants`, `initializer`, `timestepper`. When a `vegetation` component is present the
defaults switch to Richards soil hydrology and the PALADYN canopy surface stack
(`land_model.jl:117-131`). The per-step coupling order is
`atmosphere → soil → vegetation → surface_hydrology → surface_energy_balance` for auxiliaries and
`surface_hydrology → soil → vegetation` for tendencies (`land_model.jl:86-102`).

Validated on CPU by `docs/dev/2026-07/spike_land_column_cpu.jl` (soil + SEB + surface hydrology +
prescribed atmosphere + vegetation carbon; 20 steps at Δt = 60 s).

## Redundancy map: legacy physical soil/surface/climate → Terrarium config

| Legacy Agrocosm file | Terrarium replacement (constructor) | Notes / crop-specific physics to preserve |
|---|---|---|
| `parameters` `soildepth`/`layerbound` (5-layer, mm) | `ColumnGrid(arch, NF, PrescribedSpacing(Δz=...))` or `ExponentialSpacing` | Legacy layers 200/300/500/1000/1000 mm → prescribed Δz in metres. Terrarium columns are continuous-depth, not fixed 5-layer. |
| `SoilParams` texture lookup (`sand`,`silt`,`clay`,`w_sat`,`tdiff_*`; 14 soil types) | `SoilTexture(NF; sand, clay, silt)` inside a `SoilStratigraphy`/`ConstantSoilHorizon`; `SoilPorositySURFEX` + `SoilHydraulicsSURFEX` so texture drives porosity, field capacity, wilting point | Texture→quartz (`q=sand`) sets thermal conductivity; `clay` sets wilting/field capacity (SURFEX). Legacy `w_sat`/`tdiff` tables become derived, not prescribed. |
| `processes/soil/pedotransfer.jl` | `SoilHydraulicsSURFEX(NF)` (texture→hydraulics) / `ConstantSoilHydraulics` (VanGenuchten/BrooksCorey) | Legacy organic-matter effect on hydraulics → choose porosity/hydraulics scheme; organic fraction via `MineralOrganic.organic`. |
| `processes/soil/soil_temp.jl`, `water_ice_pools.jl` | `SoilThermodynamics(NF)` (enthalpy heat conduction + `FreeWater`/SFCC freeze curve) inside `SoilEnergyWaterCarbon` | Exposes auxiliary `temperature` (°C) = the crop-relevant soil temperature. Legacy 3-reservoir (wilting/available/free) partition → Terrarium freeze curve + hydraulics. |
| `processes/soil/soil_water.jl`, `infil_perc.jl` | `SoilHydrology(NF, RichardsEq(); hydraulic_properties)` | **Nitrate advective transport** currently embedded in `infil_perc.jl` is biogeochem physics — port with the soil C–N process (Phase 3), coupled to Richards water fluxes, not lost here. |
| `processes/soil/evaporation.jl` | `BareGroundEvaporation(NF; ground_resistance=SoilMoistureResistanceFactor(NF))` (bare) or `PALADYNCanopyEvapotranspiration(NF)` (canopy) | Canopy ET is the vegetated default; bare-ground evaporation is the fallback with no vegetation. |
| `processes/crop/radiation.jl` (`petpar!` daylength/PAR/eq-ET) | `PrescribedAtmosphere` inputs `surface_shortwave_down`, `daytime_length`; `DiagnosedRadiativeFluxes` | **Crop PAR / daylength astronomy** feeding photosynthesis is crop physics — port into the C3/C4 photosynthesis process (Phase 3); Terrarium supplies incoming shortwave + daytime length as inputs. |
| `processes/crop/albedo.jl` | `ConstantAlbedo(NF)` / `PrescribedAlbedo(NF)` in `SurfaceEnergyBalance` | **Canopy leaf/litter/soil/snow albedo weighting** is crop-canopy physics — if needed, add a crop `AbstractAlbedo` producing `albedo` from LAI (Phase 3); otherwise prescribe. |
| `processes/crop/interception.jl` | `PALADYNCanopyInterception(NF)` (prognostic `canopy_water`, inputs `leaf_area_index`,`stem_area_index`) | LAI/SAI must be provided by the crop vegetation process. |
| `processes/crop/transpiration.jl` | `PALADYNCanopyEvapotranspiration(NF)` + `FieldCapacityLimitedPAW(NF)` + `StaticExponentialRootDistribution` | **λ water-coupling** and demand/supply logic are crop physics coupled through `canopy_water_conductance` (stomatal) and `plant_available_water`/`soil_moisture_limiting_factor` (β). Port with the crop stomatal/photosynthesis processes (Phase 3). |
| `processes/climate/readclimate.jl` | `PrescribedAtmosphere(NF)` + `InputSources`/`InputSource` (incl. `TerrariumRastersExt`) | Inputs: `air_temperature`, `air_pressure`, `windspeed`, `specific_humidity`, `rainfall`, `snowfall`, `surface_shortwave_down`, `surface_longwave_down`, `daytime_length`, `CO2`. |
| `processes/climate/snow.jl`, `SnowParams` | **No Terrarium equivalent yet** (single/multi-layer snow) | **Genuinely missing** — flag for upstreaming (see below). |

## Crop coupling seams into Terrarium's vegetation stack

`VegetationCarbon(NF; ...)` (`Terrarium/src/processes/vegetation/vegetation_carbon.jl:6`) is a bundle
of swappable `AbstractProcess` components. Agrocosm's crop physiology replaces specific slots:

- **Photosynthesis** → subtype `AbstractPhotosynthesis{NF}`, implement `variables`,
  `compute_photosynthesis(i,j,grid,fields,photo,atmos)` (returns `(Rd, An, GPP)`) and
  `compute_auxiliary!`. Store the fields downstream consumers read: `leaf_respiration`,
  `net_assimilation`, `gross_primary_production` (`lue_photosynthesis.jl:70,412`). This is the seam
  for C3/C4 mechanistic photosynthesis + the λ solver. Inputs available: atmosphere
  `air_temperature`/`air_pressure`/`surface_shortwave_down`/`CO2`, plus
  `soil_moisture_limiting_factor`, `leaf_area_index`, `leaf_to_air_co2_ratio`.
- **Stomatal conductance** → `AbstractStomatalConductance`; produces `canopy_water_conductance`,
  `leaf_to_air_co2_ratio` (Medlyn reference at `medlyn_stomatal_conductance.jl:28,45`). The crop λ
  coupling maps here.
- **Autotrophic respiration** → `AbstractAutotrophicRespiration`; produces `autotrophic_respiration`,
  `net_primary_production` (`autotrophic_respiration.jl:32`). Crop maintenance/growth respiration.
- **Phenology** → `AbstractPhenology`; produces `phenology_factor`, `leaf_area_index`
  (`phenology.jl:22`). Crop PHU/vernalization phenology (incl. legacy `temp_stress`, `climbuf`
  vernalization) goes here.
- **Carbon dynamics / allocation** → `AbstractVegetationCarbonDynamics`; prognostic
  `carbon_vegetation`, `balanced_leaf_area_index` (`carbon_dynamics.jl:48-49`). Crop organ allocation
  + harvest index.
- **Root distribution** → `AbstractRootDistribution`; `root_density(z)` → `root_fraction`
  (`root_distribution.jl:37`). Legacy `root_distribution(beta_root)` (exponential β-profile) ports
  directly here.
- **Plant available water** → `AbstractPlantAvailableWater`; `plant_available_water`,
  `soil_moisture_limiting_factor` from texture-dependent `field_capacity`/`wilting_point`
  (`plant_available_water.jl:22,72`).

Soil C–N biogeochemistry replaces the `biogeochem` slot of `SoilEnergyWaterCarbon` (currently
`ConstantSoilCarbonDensity`) with a new `AbstractSoilBiogeochemistry` (Phase 3), coupled to soil
`temperature`, moisture, and the Richards water flux for nitrate transport.

## Genuinely-missing pieces to contribute upstream (or add in Agrocosm)

1. **Snow scheme.** Terrarium currently has no snow component in the coupled `LandModel` path; the
   soil-top heat flux is wired directly to the surface energy balance (`land_model.jl:61-64`, with an
   explicit TODO for when a snow component is added). Agrocosm's single-layer `snow.jl`/`SnowParams`
   (snowpack, insulation, cover, water-to-snow conversion) is a candidate for a Terrarium
   `AbstractSnowModel`. Multi-layer snow is on Terrarium's roadmap — upstream target.
2. **Nitrate advective transport** coupled to soil water flux (from `infil_perc.jl`) — implement as
   part of the crop soil C–N biogeochemistry process (Phase 3).
3. **Crop canopy PAR / daylength** and **LAI-dependent canopy albedo** — implement inside the crop
   photosynthesis / a crop albedo process (Phase 3) rather than upstreaming.

## Recommended crop-context soil configuration (starting point for Phase 5)

```julia
# Texture-driven soil so porosity, field capacity, and wilting point respond to sand/clay.
strat = HomogeneousSoilStratigraphy(NF;
    texture  = SoilTexture(NF; sand = 0.4, silt = 0.4, clay = 0.2),  # per-cell via InputSource in global runs
    porosity = SoilPorositySURFEX(NF))
hydrology = SoilHydrology(NF, RichardsEq(); hydraulic_properties = SoilHydraulicsSURFEX(NF))
soil = SoilEnergyWaterCarbon(NF; strat, hydrology)   # energy = SoilThermodynamics, biogeochem → crop C–N in Phase 3
# vegetation = VegetationCarbon(NF; photosynthesis = CropC3C4Photosynthesis(NF), ...)  # Phase 3
# land = LandModel(grid; soil, vegetation)
```

For gridded runs, feed per-cell `sand_fraction`/`silt_fraction`/`clay_fraction` as namespaced
`InputSource`s keyed by horizon name (see `Terrarium/test/inputs/namespaced_inputs.jl:84`), or use
`SoilGridsStratigraphy(NF)`.
