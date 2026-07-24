module Agrocosm

# Agrocosm is a downstream package built on the Terrarium.jl land-modeling
# framework. Terrarium supplies all infrastructure (grids, state, timestepping,
# architectures, parameters, I/O, diagnostics, checkpointing) and the physical
# soil/surface processes; Agrocosm contributes crop-specific processes and a
# managed-crop land model.
#
# See AGENTS.md for conventions and docs/dev/ for design notes.
using Terrarium

# Documentation and parameter macros.
using DocStringExtensions
# `@param`/`@component` are field markers consumed by `@parameterized`, not standalone macros.
using SpeedyWeatherInternals.ParameterEditing: ParameterEditing, @parameterized,
    Positive, Nonnegative, Unbounded, UnitInterval

# Reading legacy input-data files (climate/initial-condition JLD2).
import JLD2

# Framework internals used when authoring processes and kernels.
using Oceananigans.Fields: FunctionField
using Oceananigans.Utils: launch!
using Oceananigans.Operators: Δzᵃᵃᶜ
using KernelAbstractions: @kernel, @index
# Oceananigans simulation callbacks — the mechanism for the crop-management events:
# `SpecifiedTimes` for the discrete sowing/harvest jumps, `IterationInterval` for the continuous
# fertilizer application flux.
using Oceananigans: add_callback!, SpecifiedTimes, IterationInterval

# ---------------------------------------------------------------------------
# Crop parameter sets and the 12-CFT registry: infrastructure-free physics
# constants with LPJmL-derived defaults.
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
# Input data handling: build Terrarium time-varying input sources from tabulated
# climate series (reuses Oceananigans `FieldTimeSeries` + Terrarium input sources).
# ---------------------------------------------------------------------------
include("input_output/climate_forcing.jl")
export surface_climate_inputs
include("input_output/crop_initial_conditions.jl")
export load_crop_initial_conditions

# ---------------------------------------------------------------------------
# Crop processes (Terrarium-native, continuous-time), ported from the legacy
# discrete-daily LPJmL-derived implementation (now in the git history).
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
# Crop vegetation model assembling the crop processes.
# ---------------------------------------------------------------------------
include("crop/vegetation.jl")
export CropVegetation

# CFT presets: build the crop processes from the 12-CFT registry (must follow the process includes).
include("crop/cft_presets.jl")

# Top-level managed-crop model constructor (must follow CropVegetation + the soil biogeochemistry).
include("crop/crop_model.jl")
export CropModel

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

include("crop/soil_biogeochemistry.jl")
export CropSoilBiogeochemistry, soil_carbon_tendencies, soil_nitrogen_tendencies

# ---------------------------------------------------------------------------
# Crop management. Discrete lifecycle events (sowing, harvest) as documented
# Oceananigans callbacks; fertilizer as a continuous input flux.
# ---------------------------------------------------------------------------
include("crop/management.jl")
export CropCalendar, sow!, harvest!, add_crop_management!
export CropFertilization, fertilize!, add_crop_fertilization!

# The standalone LPJmL-derived reference implementation (the discrete-daily crop/soil physics and its
# bespoke infrastructure under the former `src/processes/` and `src/utils/`) was removed once its
# physics had been re-expressed as the continuous-time Terrarium processes above; it remains available
# in the git history.

end
