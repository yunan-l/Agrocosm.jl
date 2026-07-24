module Agrocosm

# Agrocosm is a downstream package built on the Terrarium.jl land-modeling
# framework. Terrarium supplies all infrastructure (grids, state, timestepping,
# architectures, parameters, I/O, diagnostics, checkpointing) and the physical
# soil/surface processes; Agrocosm contributes crop-specific processes and a
# managed-crop land model.
#
# See docs/dev/2026-07/2026-07-23_PLAN_terrarium_migration.md for the migration
# plan and AGENTS.md for conventions.
using Terrarium

# Documentation and parameter macros.
using DocStringExtensions
# `@param`/`@component` are field markers consumed by `@parameterized`, not standalone macros.
using SpeedyWeatherInternals.ParameterEditing: ParameterEditing, @parameterized,
    Positive, Nonnegative, Unbounded, UnitInterval

# Framework internals used when authoring processes and kernels.
using Oceananigans.Fields: FunctionField
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

# ---------------------------------------------------------------------------
# Crop parameter sets and the 12-CFT registry (infrastructure-free physics
# constants). These retain their LPJmL-derived defaults during the migration;
# their conversion to ModelParameters `@parameterized`/`@param` calibration
# hooks happens alongside the processes that consume them (plan Phases 3 & 5).
# ---------------------------------------------------------------------------
include("parameters/default_params.jl")
include("parameters/pft.jl")

# Numerical primitives with no framework dependencies.
include("numerics/lpj_bisect.jl")

# Global LPJmL-derived process parameters and precision handling.
export LPJmLParams, PhotoParams, ModelParameters
export SoilParams, SoilDecompParams, SoilThermalParams, SnowParams
export lpjmlparams, photoparams, soilparams, soil_decomp_params, soil_thermal_params, snowparams
export convert_precision
export FERTILIZER_MODES, fertilizer_mode

# 12-CFT crop parameter registry (LPJmL 6.1.1 order).
export PftParameters
export CROP_PFT_NAMES, CROP_PFTS, crop_pft
export cft1, cft2, cft3, cft4, cft5, cft6, cft7, cft8, cft9, cft10, cft11, cft12

# Numerics.
export lpj_bisect

# ---------------------------------------------------------------------------
# Crop processes (Terrarium-native, continuous-time). Ported incrementally in
# Phase 3 from the legacy discrete-daily implementations under src/processes/.
# ---------------------------------------------------------------------------
include("crop/root_distribution.jl")
export CropRootDistribution

include("crop/photosynthesis.jl")
export CropPhotosynthesis, C3Pathway, C4Pathway

include("crop/stomatal_conductance.jl")
export CropStomatalConductance

include("crop/carbon_dynamics.jl")
export CropCarbonDynamics

include("crop/phenology.jl")
export CropPhenology

include("crop/phenology_dynamics.jl")
export CropPhenologyDynamics, heat_unit_fraction, heat_unit_rate

include("crop/nitrogen_limitation.jl")
export CropNitrogenVcmaxLimit

include("crop/nitrogen_demand.jl")
export CropNitrogenDemand

include("crop/plant_available_water.jl")
export soil_moisture_limiting_factor, plant_available_water

include("crop/nitrogen_uptake.jl")
export CropNitrogenUptakeKinetics, nitrogen_uptake_temperature_response, root_nitrogen_uptake_potential

include("crop/soil_decomposition_response.jl")
export CropSoilDecompositionResponse, soil_decomposition_response
export soil_decomposition_temperature_response, soil_decomposition_moisture_response

include("crop/growth_respiration.jl")
export CropGrowthRespiration, growth_respiration, net_primary_production

include("crop/maintenance_respiration.jl")
export CropMaintenanceRespiration, maintenance_respiration
export maintenance_temperature_response, organ_maintenance_respiration

include("crop/harvest_index.jl")
export CropHarvestIndex, crop_harvest_index

include("crop/carbon_allocation.jl")
export CropCarbonAllocation, root_allocation_fraction, leaf_carbon_from_lai

include("crop/carbon.jl")
export CropCarbon, crop_carbon_budget

include("crop/nitrogen_allocation.jl")
export CropNitrogenAllocation, allocate_crop_nitrogen

include("crop/nitrogen.jl")
export CropNitrogen, leaf_nitrogen_limitation

# ---------------------------------------------------------------------------
# Phase 5 — crop vegetation model assembling the crop processes.
# ---------------------------------------------------------------------------
include("crop/vegetation.jl")
export CropVegetation

# CFT presets: build the crop processes from the 12-CFT registry (must follow the process includes).
include("crop/cft_presets.jl")

include("crop/soil_carbon.jl")
export CropSoilCarbon, decomposed_carbon, route_litter_carbon, heterotrophic_respiration

include("crop/nitrification.jl")
export CropNitrification, gross_nitrification
export nitrification_moisture_factor, nitrification_temperature_factor, nitrification_ph_factor

include("crop/denitrification.jl")
export CropDenitrification, gross_denitrification
export denitrification_temperature_factor, denitrification_moisture_factor

include("crop/volatilization.jl")
export CropVolatilization, ammonia_volatilization

include("crop/mineralization.jl")
export CropNitrogenMineralization, immobilization_demand, immobilization_limitation, immobilized_nitrogen

# ---------------------------------------------------------------------------
# PHASE 3+ TODO — crop and soil-biogeochemistry physics not yet ported.
#
# The following source files remain on disk as the reference implementation but
# are NOT included yet: they depend on the deleted standalone infrastructure
# (`launch_1D!`, `Parameters.@unpack`, the legacy state containers) and must be
# re-expressed as continuous-time Terrarium `AbstractProcess` implementations
# with `variables()`/`compute_auxiliary!`/`compute_tendencies!`.
#
#   Crop physiology (Phase 3):
#     processes/crop/{photosynthesis,lambda_solver,respiration,carbon_allocation,
#                     crop_carbon,lai_crop,phenology}.jl
#     processes/crop/{nitrogen_allocation,nitrogen_demand,nitrogen_uptake,
#                     nitrogen_vcmax_limit}.jl
#     processes/climate/{temp_stress,climbuf,spinup_climbuf}.jl  (temp stress +
#                     vernalization are crop physics, not framework climate)
#   Surface coupling to reuse Terrarium processes (Phase 2/3):
#     processes/crop/{radiation,albedo,interception,transpiration}.jl
#   Soil C–N biogeochemistry (Phase 3):
#     processes/soil/{soil_carbon,soil_nitrogen,nitrogen_transform,soil_response,
#                     litter_routing,surface_litter}.jl
#   Physical soil to replace with Terrarium soil processes (Phase 2):
#     processes/soil/{soil_temp,water_ice_pools,soil_water,infil_perc,
#                     pedotransfer,evaporation}.jl  and  processes/climate/{readclimate,snow}.jl
#   Management events (Phase 4, documented discrete-time exceptions):
#     processes/crop/{cultivate,harvesting,fertilizer}.jl, processes/soil/tillage.jl
#     utils/tools.jl (fixed-365 calendar convention)
#   State schemas — reference for Terrarium `variables()` (Phase 3):
#     processes/initialization/**
# ---------------------------------------------------------------------------

end
