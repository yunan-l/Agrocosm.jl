module Agrocosm

# Write your package code here.
# NUMERICS
using Statistics

# GPU PARALLEL
# import KernelAbstractions: @kernel, @index, @inbounds # get_backend, synchronize
using KernelAbstractions # GPU/CPU parallelization
using CUDA, Adapt

# INPUT OUTPUT
using DataFrames, NCDatasets, Random, Dates
import JLD2: @load, @save, load

# PARAMETER HANDLING
import Parameters: @with_kw, @unpack
import MuladdMacro: @muladd

# STRUCTURES
export LPJmLParams, PftParameters, PhotoParams, ModelParameters, PetPar, Output
export DailyWeather, ClimBuf, CO2
export Crop, CropState, CropFluxes, CropAuxiliary, CropWorkspace
export CropPhenology, CropCanopyState, CropCarbonState, CropNitrogenState, CropWaterState
export CropCalendarState, CropCanopyAuxiliary, CropPhotosynthesisAuxiliary
export CropCarbonFluxes, CropNitrogenFluxes, CropWaterFluxes, CropEvents
export CropStressAuxiliary, ManagedLand
export SoilParams, SoilDecompParams, SoilThermalParams, SnowParams, Soil
export SoilProperties, SoilWater, SoilThermal, SoilCarbon, SoilNitrogen
export SoilDecomposition, SoilManagement, SoilSurfaceLitter, SoilSnow
export CropOutput, SoilOutput, ClimateOutput, CalendarOutput

# PARAMETERS (PFTs)
export lpjmlparams, photoparams, soilparams, soil_decomp_params, soil_thermal_params, snowparams, cft1, cft2, cft3, cft4
export convert_precision

# INITIALIZATION
export init_states!, init_climbuf, init_crop, init_pet, init_soil, init_output
export init_weather, init_managed_land
export WaterBalance, init_water_balance
export NitrogenBalance, init_nitrogen_balance
export CarbonBalance, init_carbon_balance
export ThermalBalance, init_thermal_balance

# CLIMATE
export annual_climbuf!, daily_climbuf!, infil_perc!, spin_up_climbuf!, update_climbuf!, readclimate!, snow!

# PHYSICS FUNCTIONS
# RADIATION
export albedo!, petpar!, apar_crop!, apar_crop_maize!

# CROP
export photosynthesis_C3!, photosynthesis_C4!, carbon_allocation!, respiration!
export phenology_crop!, lai_crop!, lai_deficit!, cultivate!, harvest_crop!, fertilizer!
export transpiration!, interception!
export crop_nitrogen!, ndemand_crop!, nuptake_crop!
export limit_vmax_by_nitrogen!
export root_distribution, temp_stress
export lpj_bisect, solve_lambda_c3_lpj, solve_lambda_c4_lpj
export solve_lambda_c3!, solve_lambda_c4!
export crop_carbon!, crop_carbon_hybrid!, hybrid_photos_C3!, hybrid_photos_C4!
export waterlogging_stress!

# SOIL
export apply_percolation_enthalpy!, soil_temperature!
export pedotransfer!, soil_carbon!
export evaporation!, soil_infiltration!, soil_evapotranspiration!
export soil_nitrogen!, nitrogen_transform!, soil_cn_decomposition!, post_crop_nitrogen_losses!
export soil_decomp_response!
export update_surface_litter_properties!, surface_litter_interception!
export litter_tillage!, litter_bioturbation!

# UNITS
export deg2rad, ppm2Pa, ppm2bar, hour2day, hour2sec, degCtoK

# DATA
export InitialDataLoader, ClimateDataLoader, DataLoader, DataLoader_winter_wheat
export write_output_nc

# DAILY CROP SIMULATIONS
export daily_crop_C3!, daily_crop_C4!
export CropSimulation, initialize_simulation, run_simulation!, simulation_summary


# process-based crop model
# Parameters
include("parameters/default_params.jl")
include("parameters/pft.jl")

# Numerics
include("numerics/lpj_bisect.jl")

# Initialization
include("processes/initialization/climate/climate.jl")
include("processes/initialization/management/managed_land.jl")
include("processes/initialization/output/output.jl")
include("processes/initialization/crop/phenology.jl")
include("processes/initialization/crop/canopy.jl")
include("processes/initialization/crop/carbon.jl")
include("processes/initialization/crop/nitrogen.jl")
include("processes/initialization/crop/water.jl")
include("processes/initialization/crop/calendar.jl")
include("processes/initialization/crop/photosynthesis.jl")
include("processes/initialization/crop/crop.jl")
include("processes/initialization/soil/properties.jl")
include("processes/initialization/soil/water.jl")
include("processes/initialization/soil/thermal.jl")
include("processes/initialization/soil/carbon.jl")
include("processes/initialization/soil/nitrogen.jl")
include("processes/initialization/soil/decomposition.jl")
include("processes/initialization/soil/management.jl")
include("processes/initialization/soil/surface_litter.jl")
include("processes/initialization/soil/snow.jl")
include("processes/initialization/soil/soil.jl")
include("processes/initialization/init_states.jl")

# Diagnostics
include("diagnostics/water_balance.jl")
include("diagnostics/nitrogen_balance.jl")
include("diagnostics/carbon_balance.jl")
include("diagnostics/thermal_balance.jl")

# Climate
include("processes/climate/climbuf.jl")
include("processes/climate/temp_stress.jl")
include("processes/climate/spinup_climbuf.jl")
include("processes/climate/readclimate.jl")
include("processes/climate/snow.jl")

# Crop
include("processes/crop/cultivate.jl")
include("processes/crop/phenology.jl")
include("processes/crop/photosynthesis.jl")
include("processes/crop/lambda_solver.jl")
include("processes/crop/carbon_allocation.jl")
include("processes/crop/crop_carbon.jl")
include("processes/crop/lai_crop.jl")
include("processes/crop/radiation.jl")
include("processes/crop/albedo.jl")
include("processes/crop/respiration.jl")
include("processes/crop/interception.jl")
include("processes/crop/transpiration.jl")
include("processes/crop/nitrogen_allocation.jl")
include("processes/crop/nitrogen_demand.jl")
include("processes/crop/nitrogen_uptake.jl")
include("processes/crop/nitrogen_vmax_limit.jl")
include("processes/crop/fertilizer.jl")
include("processes/crop/harvesting.jl")
include("processes/crop/waterlogging_stress.jl")

# Soil
include("processes/soil/water_ice_pools.jl")
include("processes/soil/pedotransfer.jl")
include("processes/soil/evaporation.jl")
include("processes/soil/soil_temp.jl")
include("processes/soil/surface_litter.jl")
include("processes/soil/litter_routing.jl")
include("processes/soil/nitrogen_transform.jl")
include("processes/soil/infil_perc.jl")
include("processes/soil/soil_water.jl")
include("processes/soil/soil_carbon.jl")
include("processes/soil/soil_nitrogen.jl")
include("processes/soil/soil_response.jl")

# Input and output
include("input_output/climate_data_loader.jl")
include("input_output/initial_data_loader.jl")
include("input_output/write_output_nc.jl")

# Utilities
include("utils/kernel_launch.jl")
include("utils/visualization.jl")
include("utils/conversions.jl")
include("utils/load_nc.jl")
include("utils/tools.jl")

# Daily crop simulations
include("simulations/daily_crop_C3.jl")
include("simulations/daily_crop_C4.jl")
include("simulations/simulation_api.jl")

end
