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
import JLD2: @load, @save

# PARAMETER HANDLING
import Parameters: @with_kw, @unpack
import MuladdMacro: @muladd

# STRUCTURES
export LPJmLParams, PftParameters, PhotoParams, Photos, PetPar, Output
export DailyWeather, ClimBuf, CO2
export Crop, Calendar, Managed_land
export SoilParams, SoilDecompParams, SnowParams, Soil

# PARAMETERS (PFTs)
export lpjmlparams, photoparams, soilparams, soil_decomp_params, snowparams, cft1, cft2, cft3, cft4

# INITIALIZATION
export init_states!, init_climbuf, init_crop, init_pet, init_soil, init_output
export WaterBalance, init_water_balance

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
export root_distribution, temp_stress
export crop_carbon!, crop_carbon_hybrid!, hybrid_photos_C3!, hybrid_photos_C4!
export waterlogging_stress!

# SOIL
export soiltemp_lag!
export pedotransfer!, soil_carbon!, update_litc_tillage!, update_lit_winter_wheat!
export evaporation!, soil_water!
export soil_nitrogen!, nitrogen_transform!, update_litn_tillage!
export soil_decomp_response!

# UNITS
export deg2rad, ppm2Pa, ppm2bar, hour2day, hour2sec, degCtoK

# DATA
export InitialDataLoader, ClimateDataLoader, DataLoader, DataLoader_winter_wheat
export write_output_nc

# DAILY CROP SIMULATIONS
export daily_crop_C3!, daily_crop_C4!


# process-based crop model
# Parameters
include("parameters/default_params.jl")
include("parameters/pft.jl")

# Initialization
include("processes/initialization/define_structs.jl")
include("processes/initialization/init_states.jl")
include("processes/initialization/init_structs.jl")

# Diagnostics
include("diagnostics/water_balance.jl")

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
include("processes/crop/fertilizer.jl")
include("processes/crop/harvesting.jl")
include("processes/crop/waterlogging_stress.jl")

# Soil
include("processes/soil/pedotransfer.jl")
include("processes/soil/evaporation.jl")
include("processes/soil/soil_temp.jl")
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
include("utils/callback.jl")
include("utils/tools.jl")

# Daily crop simulations
include("simulations/daily_crop_C3.jl")
include("simulations/daily_crop_C4.jl")

end
