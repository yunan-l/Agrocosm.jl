# Agrocosm variable inventory

This document is the authoritative lifecycle inventory for model variables. It
complements the LPJmL-style field comments beside every Julia struct field.

## Classification rules

| Class | Definition | Daily lifecycle | Restart/output rule |
|---|---|---|---|
| `state` | Prognostic quantity whose absence changes a resumed model transition | Read from the previous day and updated | Must be checkpointed |
| `flux` | Transfer occurring during the current daily time step | Must be overwritten on every relevant path | Not checkpointed; may be output or used in balance ledgers |
| `auxiliary` | Algebraic diagnostic, forcing-derived quantity, static spatial coefficient, or explicit process memory | Recomputed, overwritten, or deliberately retained according to its field description | Recomputable values are excluded from restart; retained process memory is listed explicitly |
| `event` | One-day discrete transition | Zero except on its event day | Not checkpointed; may be output |
| `workspace` | Preallocated implementation scratch with no scientific meaning | Arbitrary after a kernel completes | Never checkpointed or exposed as scientific output |
| `forcing/configuration` | External driver, management input, parameter, geometry, precision, or backend choice | Supplied or buffered rather than predicted | Re-read or stored as simulation configuration |
| `output` | Time-series copy selected for users | Appended after completed steps | Scientific result, not restart state |
| `output bookkeeping` | In-progress output aggregate that never feeds a process | Updated until its reporting boundary | Save with a full simulation checkpoint to preserve incomplete reports |
| `balance diagnostic` | Optional conservation ledger | Recorded for each completed day | Debug/validation result, not restart state |

Array conventions are `(cell)`, `(layer, cell)`, or `(day, cell)`. The three
litter rows are surface litter, tillage-incorporated litter, and root litter.
The default five soil layers have thicknesses 200, 300, 500, 1000, and 1000 mm.

## Crop state

These fields form `crop.state` and are the prognostic crop checkpoint.

| Path | Meaning | Unit |
|---|---|---|
| `crop.state.phenology.vdsum` | Accumulated effective vernalization | day equivalent |
| `crop.state.phenology.husum` | Accumulated heat units since sowing | °C day |
| `crop.state.phenology.senescence` | Current senescence mode | Bool |
| `crop.state.phenology.senescence_previous` | Previous-day senescence mode used at transitions | Bool |
| `crop.state.phenology.harvesting` | Current harvest-readiness mode | Bool |
| `crop.state.phenology.harvesting_previous` | Previous-day harvest-readiness mode | Bool |
| `crop.state.phenology.growing_days` | Days since cultivation | day |
| `crop.state.phenology.is_growing` | Active crop-stand mode | 0/1 |
| `crop.state.canopy.lai` | Potential phenological leaf-area index | m² leaf m⁻² ground |
| `crop.state.canopy.laimax_adjusted` | LAI retained at senescence onset | m² leaf m⁻² ground |
| `crop.state.canopy.lai_npp_deficit` | LAI unsupported by available carbon | m² leaf m⁻² ground |
| `crop.state.carbon.biomass` | Total live crop carbon | gC m⁻² |
| `crop.state.carbon.leaf` | Live leaf carbon | gC m⁻² |
| `crop.state.carbon.root` | Live root carbon | gC m⁻² |
| `crop.state.carbon.pool` | Mobile/intermediate carbon pool | gC m⁻² |
| `crop.state.carbon.storage` | Harvestable storage-organ carbon | gC m⁻² |
| `crop.state.nitrogen.total` | Total live crop nitrogen | gN m⁻² |
| `crop.state.nitrogen.leaf` | Leaf nitrogen | gN m⁻² |
| `crop.state.nitrogen.root` | Root nitrogen | gN m⁻² |
| `crop.state.nitrogen.pool` | Mobile/intermediate nitrogen pool | gN m⁻² |
| `crop.state.nitrogen.storage` | Storage-organ nitrogen | gN m⁻² |
| `crop.state.nitrogen.pending_manure` | Manure N awaiting split application | gN m⁻² |
| `crop.state.nitrogen.pending_fertilizer` | Mineral fertilizer N awaiting split application | gN m⁻² |
| `crop.state.nitrogen.stress_sum` | Seasonal integral of N sufficiency | dimensionless day |
| `crop.state.nitrogen.sufficiency` | Prior-day N sufficiency used by next-day LAI growth | 0–1 |
| `crop.state.water.demand_sum` | Seasonal transpiration-demand integral | mm |
| `crop.state.water.supply_sum` | Seasonal transpiration-supply integral | mm |
| `crop.state.water.sufficiency` | Prior-day water sufficiency used by next-day LAI growth | 0–1 |

`biomass` is retained for LPJmL-compatible process flow, although it should be
consistent with the organ-carbon pools and is checked as a derived total.

### Why each crop group is, or is not, prognostic

| Variables | Classification | Checkpoint criterion |
|---|---|---|
| `vdsum`, `husum`, phase flags, `growing_days`, `is_growing` | prognostic state | The next phenology transition needs their previous values and cannot reconstruct the season phase from forcing alone. |
| `lai`, `laimax_adjusted`, `lai_npp_deficit` | prognostic state | Before daily phenology runs, interception/radiation use prior canopy size; senescence and carbon-deficit history cannot be recovered from forcing alone. |
| carbon and nitrogen organ pools, `biomass` | prognostic state | These are living stocks advanced by daily C/N transfers. |
| pending manure/fertilizer and `stress_sum` | prognostic state | They encode incomplete split application and seasonal stress history. |
| N/water `sufficiency`, `demand_sum`, `supply_sum` | prognostic state | The next LAI and allocation update use their previous values/integrals. |
| `phu`, `winter_type`, `sowing_date`, root distribution | configuration/process memory | They are supplied by initialization and do not evolve as a crop response. |
| `fphu`, `phenology_fraction`, actual LAI, radiation, conductance, stress diagnostics, photosynthesis variables | diagnostic | They are reconstructed within the daily process chain from state, configuration, and forcing. `phenology_fraction` is a kernel-local scalar, not an allocated field. |
| C/N/water fluxes | flux | They represent only the current transition and owner kernels overwrite them. |
| sowing/harvest masks | event | They are one-day discrete transitions. |
| `output.annual.yield`, `output.annual.harvest_date` | output bookkeeping | They preserve an incomplete annual report but never alter a crop or soil process. |

## Crop daily fluxes

| Path | Meaning | Unit |
|---|---|---|
| `crop.fluxes.carbon.yield` | Storage-organ carbon harvested today | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.harvest_export` | Total crop carbon exported at harvest | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.npp` | GPP minus leaf, root, organ-maintenance, and growth respiration | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.respiration` | Total plant respiration | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.gross_assimilation` | Gross canopy assimilation/GPP | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.net_assimilation` | Nonnegative daytime assimilation after leaf respiration | gC m⁻² day⁻¹ |
| `crop.fluxes.carbon.water_limited_assimilation` | Assimilation expressed as transpiration demand | mm day⁻¹ |
| `crop.fluxes.carbon.leaf_respiration` | Leaf/dark respiration | gC m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.uptake` | Total mineral N transferred to crop | gN m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.auto_fertilizer` | Uptake supplied by automatic fertilizer | gN m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.seed_input` | Nitrogen introduced with seed | gN m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.prescribed_manure_input` | Scheduled manure N applied today | gN m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.prescribed_fertilizer_input` | Scheduled mineral N applied today | gN m⁻² day⁻¹ |
| `crop.fluxes.nitrogen.harvest_export` | Plant N removed at harvest | gN m⁻² day⁻¹ |
| `crop.fluxes.water.interception` | Canopy-intercepted water evaporated today | mm day⁻¹ |
| `crop.fluxes.water.transpiration_layer` | Root uptake/transpiration by soil layer | mm day⁻¹ |

## Crop auxiliary variables and process memory

| Path | Meaning | Lifecycle |
|---|---|---|
| `crop.auxiliary.phenology.phu` | Heat units required for maturity (°C day) | Seasonal process memory |
| `crop.auxiliary.phenology.winter_type` | Vernalization requirement flag | Seasonal process memory |
| `crop.auxiliary.phenology.fphu` | Current heat-unit fraction reconstructed as `husum / phu` | Daily diagnostic |
| `crop.auxiliary.calendar.sowing_date` | Prescribed sowing day of year | Configuration/process memory |
| `crop.auxiliary.canopy.actual_lai` | Nonnegative carbon-supported LAI | Daily derived diagnostic |
| `crop.auxiliary.canopy.flaimax` | Current-day fraction of maximum LAI from phenology | Daily derived diagnostic |
| `crop.auxiliary.canopy.albedo` | Effective crop-covered surface albedo, including canopy, litter, soil, and snow | Daily derived diagnostic |
| `crop.auxiliary.canopy.fpar` | Fraction of PAR absorbed | Daily derived diagnostic |
| `crop.auxiliary.canopy.apar` | Absorbed PAR (J m⁻² day⁻¹) | Daily derived diagnostic |
| `crop.auxiliary.canopy.canopy_conductance` | Bulk canopy conductance (mm s⁻¹) | Daily derived diagnostic |
| `crop.auxiliary.canopy.canopy_wet` | Wet-canopy evaporation fraction | Daily derived diagnostic |
| `crop.auxiliary.photosynthesis.potential_vcmax` | Potential carboxylation capacity | Daily derived diagnostic |
| `crop.auxiliary.photosynthesis.vcmax` | Realized carboxylation capacity | Daily derived diagnostic |
| `crop.auxiliary.photosynthesis.nitrogen_limitation` | Realized/potential Vcmax ratio | Daily diagnostic |
| `crop.auxiliary.photosynthesis.lambda` | Intercellular/ambient CO₂ ratio | Daily solved diagnostic |
| `crop.auxiliary.photosynthesis.temperature_stress` | Photosynthetic temperature multiplier | Daily derived diagnostic |
| `crop.auxiliary.stress.nitrogen_demand_total` | Potential whole-plant N demand | Daily diagnostic |
| `crop.auxiliary.stress.nitrogen_demand_leaf` | Potential leaf N demand | Daily diagnostic |
| `crop.auxiliary.stress.nitrogen_deficit` | Unmet crop N demand | Daily diagnostic |
| `crop.auxiliary.stress.water_deficit` | Allocation water-deficit factor | Daily diagnostic |
| `crop.auxiliary.root.distribution` | Static fraction of roots by layer | Derived from PFT parameter `beta_root` |
| `crop.auxiliary.root.zone_available_water` | Top-three-layer root-weighted plant-available water | mm; daily diagnostic |

`crop.events.sowing` and `crop.events.harvest` are one-day 0/1 events.
`CropWorkspace` is currently empty. The scientific restart contains
`crop.state` plus static auxiliary configuration required to reconstruct the
next transition; daily fluxes, events, and recomputable auxiliaries are
excluded.

## Output bookkeeping

`output.annual.yield` and `output.annual.harvest_date` retain a partial
calendar-year report between a harvest and day 365. They never feed crop,
soil, or management processes, so they are not prognostic crop state. A full
simulation checkpoint must retain this small accumulator together with output
arrays whenever a run may resume before the year-end report is emitted.

## Soil variable classification

The soil structs remain process-grouped, so the lifecycle class is attached to
each field rather than inferred from the owning struct.

### Soil state

| Path | Meaning | Unit |
|---|---|---|
| `soil.water.storage` | Liquid water by layer | mm |
| `soil.water.ice_storage` | Total ice by layer | mm water equivalent |
| `soil.water.wilting_ice_fraction` | Frozen fraction of wilting-point water | 0–1 |
| `soil.water.available_ice_storage` | Ice in plant-available water | mm |
| `soil.water.free_ice_storage` | Ice in gravitational water | mm |
| `soil.thermal.temperature` | Layer temperature | °C |
| `soil.thermal.enthalpy` | Volumetric layer enthalpy relative to 0 °C | J m⁻³ |
| `soil.thermal.frozen_fraction` | Fraction of layer water frozen | 0–1 |
| `soil.thermal.freeze_depth` | Effective frozen depth in layer | mm |
| `soil.thermal.water_reference` | Phase-change reference water stock | mm |
| `soil.thermal.initialized` | Thermal-profile initialization mode | Bool |
| `soil.carbon.litter` | C in surface/incorporated/root litter | gC m⁻² |
| `soil.carbon.fast` | Fast SOC by layer | gC m⁻² |
| `soil.carbon.slow` | Slow SOC by layer | gC m⁻² |
| `soil.nitrogen.nitrate` | Soil NO₃-N by layer | gN m⁻² |
| `soil.nitrogen.ammonium` | Soil NH₄-N by layer | gN m⁻² |
| `soil.nitrogen.litter` | Organic N in three litter classes | gN m⁻² |
| `soil.nitrogen.fast` | Fast organic-soil N by layer | gN m⁻² |
| `soil.nitrogen.slow` | Slow organic-soil N by layer | gN m⁻² |
| `soil.management.tillage_density_factor` | Tilled-topsoil bulk density relative to settled soil | 0–1 |
| `soil.surface_litter.dry_matter` | Surface litter dry matter | gDM m⁻² |
| `soil.surface_litter.depth` | Effective litter thickness | m |
| `soil.surface_litter.cover` | Fractional litter cover | 0–1 |
| `soil.surface_litter.water_storage` | Water retained by litter | mm |
| `soil.surface_litter.temperature` | Litter temperature | °C |
| `soil.snow.pack` | Snow water-equivalent stock | mm |

### Soil daily fluxes

| Owner | Fields | Unit |
|---|---|---|
| `soil.water` | `evaporation`, `influx`, `outflux`, `surface_runoff`, `lateral_runoff`, `bottom_drainage`, `infiltration`, `percolation` | mm day⁻¹ |
| `soil.thermal` | `percolation_energy`, `surface_energy_flux`, `untracked_water_energy_flux`, `rain_energy_input`, `snowmelt_energy_input`, `lateral_runoff_energy_output`, `bottom_drainage_energy_output` | J m⁻² day⁻¹ |
| `soil.thermal` | `energy_residual`, `percolation_energy_residual` | J m⁻² |
| `soil.carbon` | `input`, `decomposed_litter`, `decomposed_fast`, `decomposed_slow`, `litter_to_fast`, `litter_to_slow`, `heterotrophic_respiration` | gC m⁻² day⁻¹ |
| `soil.nitrogen` | `input`, `decomposed_litter`, `decomposed_fast`, `decomposed_slow`, `litter_to_fast`, `litter_to_slow`, `mineralization`, `immobilization`, `nitrification`, `n2o_nitrification`, `denitrification`, `n2o_denitrification`, `n2_denitrification`, `volatilization`, `leaching` | gN m⁻² day⁻¹ |
| `soil.management` | `tillage_carbon`, `bioturbation_carbon` | gC m⁻² day⁻¹ |
| `soil.management` | `tillage_nitrogen`, `bioturbation_nitrogen` | gN m⁻² day⁻¹ |
| `soil.surface_litter` | `interception`, `evaporation` | mm day⁻¹ |
| `soil.snow` | `melt`, `sublimation`, `runoff` | mm day⁻¹ |

### Soil auxiliary, properties, and workspace

| Owner | Fields | Class/meaning |
|---|---|---|
| `soil.properties` | `sand_fraction`, `clay_fraction`, `ph`, `layer_depth` | Static forcing/configuration |
| `soil.water` | `relative_content`, `free_water` | Daily hydraulic diagnostics |
| `soil.water` | `wilting_fraction`, `wilting_storage`, `field_capacity`, `saturation_fraction`, `saturation_storage`, `beta`, `holding_capacity_fraction`, `holding_capacity_storage`, `saturated_conductivity` | Daily derived hydraulic properties; top layer responds to tillage density |
| `soil.thermal` | `heat_capacity_frozen`, `heat_capacity_unfrozen`, `latent_heat`, `conductivity_frozen`, `conductivity_unfrozen` | Daily thermal properties |
| `soil.thermal` | `diffusivity_0`, `diffusivity_15` | Static soil-type parameters |
| `soil.carbon` | `litter_response` | Daily environmental response |
| `soil.nitrogen` | `litter_response` | Daily environmental response |
| `soil.decomposition` | `response`, `litter_response` | Daily temperature/moisture response |
| `soil.decomposition` | `shift_fast`, `shift_slow` | Fixed post-spin-up C/N input-routing configuration; not checkpoint state |
| `soil.decomposition` | `layer_scratch_1`, `layer_scratch_2`, `surface_scratch_1`, `surface_scratch_2` | Workspace; excluded from restart/output |
| `soil.management.tillage_fraction` | Litter-class routing matrix | Static management coefficient |
| `soil.surface_litter` | `water_capacity`, `conductivity` | Daily derived physical properties |
| `soil.snow` | `height`, `fraction` | Derived snow-surface diagnostics |

## Climate and management forcing

| Struct | Fields |
|---|---|
| `DailyWeather` | `temp`, `prec`, `swr`, `lwr`, `wind`, `daily_co2`, `annual_co2` |
| `PetPar` | `daylength`, `par`, `eeq`, `albedo` |
| `ClimBuf` | `temp`, `mtemp`, `mtemp20`, `min_temp`, `atemp`, `atemp_mean`, `V_req_a`, `V_req` |
| `ManagedLand` | `manure`, `fertilizer`, `residue_fraction`, `latitude` |

`ClimBuf` is forcing-derived **checkpoint state**: it is not an ecosystem
stock, but its rolling histories are read by future phenology and soil-thermal
steps and therefore must be preserved by a restart. `PetPar` and
`DailyWeather` are overwritten forcing buffers.

## Output variables

`CropOutput`: `gpp`, `npp`, `lambda`, `potential_vcmax`, `vcmax`,
`nitrogen_limitation`, `respiration`, `biomass`, `lai`, `storage_carbon`,
`yield`, `vegetation_carbon`, `vegetation_nitrogen`, `fphu`, `water_deficit`,
and `growing_mask`.

`SoilOutput`: `ecosystem_respiration`, `litter_carbon`, `fast_carbon`,
`slow_carbon`, `water_storage`, `litter_nitrogen`, `fast_nitrogen`,
`slow_nitrogen`, `heterotrophic_respiration`, and `evapotranspiration`.

`ClimateOutput`: `equilibrium_evapotranspiration`, `precipitation`, and
`temperature`. `CalendarOutput`: `harvesting_mask`, `harvesting_year`,
`harvest_date`, `sowing_event`, and `harvest_event`.

Output rows correspond only to completed simulation days; initial state is not
stored as a synthetic day-zero row.

## Balance diagnostics

All diagnostic fields have shape `(day, cell)` and do not feed back into the
simulation.

| Ledger | Fields |
|---|---|
| `CarbonBalance` | `plant_before`, `plant_after`, `soil_before`, `soil_after`, `total_before`, `total_after`, `net_primary_production`, `seed_input`, `manure_input`, `residue_transfer`, `harvest_export`, `heterotrophic_respiration`, `litter_respiration`, `fast_pool_respiration`, `slow_pool_respiration`, `residual`, `relative_residual` |
| `NitrogenBalance` | `plant_before`, `plant_after`, `mineral_before`, `mineral_after`, `organic_before`, `organic_after`, `total_before`, `total_after`, `root_uptake`, `seed_input`, `prescribed_fertilizer_input`, `prescribed_manure_input`, `automatic_fertilizer_input`, `harvest_export`, `mineralization`, `immobilization`, `nitrification`, `n2o_nitrification`, `denitrification`, `n2o_denitrification`, `n2_denitrification`, `volatilization`, `gaseous_loss`, `leaching_loss`, `residual`, `relative_residual` |
| `WaterBalance` | `precipitation`, `rain_after_snow`, `soil_storage_before`, `soil_storage_after`, `soil_ice_storage_before`, `soil_ice_storage_after`, `snow_storage_before`, `snow_storage_after`, `litter_storage_before`, `litter_storage_after`, `snowmelt`, `snow_sublimation`, `snow_runoff`, `unaccounted_snow_flux`, `interception`, `litter_interception`, `litter_evaporation`, `transpiration`, `evaporation`, `surface_runoff`, `lateral_runoff`, `bottom_drainage`, `remaining_infiltration`, `residual` |
| `ThermalBalance` | `surface_energy_flux`, `energy_residual`, `untracked_water_energy_flux`, `rain_energy_input`, `snowmelt_energy_input`, `lateral_runoff_energy_output`, `bottom_drainage_energy_output`, `percolation_energy_residual`, `column_energy`, `total_ice_storage`, `wilting_ice_storage`, `available_ice_storage`, `free_ice_storage`, `ice_pool_residual`, `maximum_frozen_fraction`, `minimum_temperature`, `maximum_temperature` |

## Parameters and configuration

Parameter fields are configuration and never prognostic state. The source
declarations contain the units and detailed LPJmL-style field comments.

### Crop/PFT parameters

| Group | Fields |
|---|---|
| Identity/pathway | `name`, `plant_type`, `path` |
| Temperature envelopes | `temp.low`, `temp.high`, `temp_co2.low`, `temp_co2.high`, `temp_photos.low`, `temp_photos.high` |
| Vernalization | `tv_eff.low`, `tv_eff.high`, `tv_opt.low`, `tv_opt.high`, `pvd_max` |
| Photoperiod | `psens`, `pb`, `ps` |
| Phenology and canopy | `basetemp.low`, `basetemp.high`, `fphuc`, `flaimaxc`, `fphuk`, `flaimaxk`, `fphusen`, `flaimaxharvest`, `laimax`, `laimin`, `hlimit` |
| Carbon and radiation traits | `b`, `albedo_leaf`, `albedo_litter`, `alphaa`, `lightextcoeff`, `longevity`, `sla`, `respcoeff`, `shapesenescencenorm`, `fpc` |
| Organ stoichiometry | `nc_ratio.root`, `nc_ratio.sto`, `nc_ratio.pool`, `ratio.root`, `ratio.sto`, `ratio.pool`, `ncleaf.low`, `ncleaf.median`, `ncleaf.high` |
| Litter and roots | `k_litter10.leaf`, `k_litter10.root`, `beta_root` |
| Water exchange | `intc`, `emax`, `gmin` |
| N storage/uptake | `knstore`, `no3_uptake.vmax`, `no3_uptake.kmin`, `no3_uptake.Km`, `nh4_uptake.vmax`, `nh4_uptake.kmin`, `nh4_uptake.Km` |
| Harvest index | `hiopt`, `himin` |

### Global LPJmL-derived parameters

| Group | `LPJmLParams` fields |
|---|---|
| Photosynthesis | `ko25`, `kc25`, `theta`, `alphac3`, `alphac4`, `k`, `LAMBDA_OPT` |
| Plant respiration | `r_growth` |
| Soil decomposition | `e0`, `temp_response`, `k_soil10.fast`, `k_soil10.slow`, `fastfrac`, `atmfrac` |
| Litter routing | `residue_frac`, `bioturbate` |
| Evaporation/hydrology | `ALPHAM`, `GM`, `PRIESTLEY_TAYLOR`, `soildepth_evap`, `p`, `soil_infil`, `soil_infil_litter`, `percthres`, `maxsnowpack` |
| Soil physics | `MINERALDENS` |
| Plant N uptake | `k_temp`, `T_0`, `T_m`, `T_r` |
| Soil N transformations | `k_max`, `k_2`, `soil_cn_ratio`, `immobilization_k`, `nitrification_a`, `nitrification_b`, `nitrification_c`, `nitrification_d`, `CDN`, `n2o_denit_frac` |
| Volatilization | `volatil_wind`, `volatil_length` |
| Fertilizer/manure | `manure_cn`, `nfert_split_frac`, `nmanure_nh4_frac`, `nfert_no3_frac` |

### Other physical and numerical parameters

| Struct | Fields |
|---|---|
| `PhotoParams` | `po2`, `p`, `q10ko`, `q10kc`, `q10tau`, `tau25`, `cmass`, `cq`, `lambdamc4`, `lambdamc3`, `tmc3`, `tmc4` |
| `SoilParams` | `sand`, `silt`, `clay`, `w_sat`, `tdiff_0`, `tdiff_15`, `soildepth` |
| `SnowParams` | `tsnow`, `snow_skin_depth`, `th_diff_snow`, `lambda_snow`, `c_water2ice`, `c_watertosnow`, `c_roughness` |
| `SoilThermalParams` | `seconds_per_day`, `diffusivity_conversion`, `soil_heat_capacity`, `litter_carbon_fraction`, `litter_bulk_density`, `litter_porosity`, `litter_conductivity_dry`, `litter_conductivity_saturated_unfrozen`, `litter_conductivity_saturated_frozen`, `mineral_heat_capacity`, `water_heat_capacity`, `ice_heat_capacity`, `volumetric_fusion_heat`, `solid_conductivity`, `water_conductivity`, `ice_conductivity`, `phase_change_substeps` |
| `SoilDecompParams` | `e0`, `intercept`, `moist3`, `moist2`, `moist1`, `eps` |

`ModelParameters` groups `LPJmLParams`, `PhotoParams`, `SnowParams`,
`SoilThermalParams`, and `SoilDecompParams` in the selected floating-point
precision. `SoilParams` is the soil-type lookup table. These values are
configuration, never prognostic state.

`CropSimulation` stores `pft`, the runtime `state` container, optional
`diagnostics`, `model_parameters`, run `config`, and `simulated_days`. Its
runtime `state` named tuple contains crop/soil state together with forcing
buffers and output for API convenience; that API grouping must not be confused
with the scientific lifecycle classification in this document.
