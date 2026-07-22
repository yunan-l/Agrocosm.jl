module AgrocosmData

using Dates
using NCDatasets
using TOML

include("contracts.jl")
include("catalog.jl")
include("grid.jl")
include("netcdf.jl")
include("masks.jl")
include("soil.jl")
include("management.jl")
include("climate.jl")

export DATA_SCHEMA_VERSION
export DataProvenance, DatasetSpec, DatasetCatalog, PFTRegistry, ManagementBands
export GridIndex, CellSelection, CompactVariable, TimeCellData, CropMask
export CO2Series, ClimateBlock, ClimateBlockReader
export SoilLookup, SoilData, DEFAULT_SOIL_LOOKUP_VERSION
export load_catalog, dataset, pft_index, pft_name
export read_grid, all_cells, select_cells, compact_spatial, expand_to_grid
export read_compact_variable, read_static_cell
export build_crop_mask
export default_soil_lookup, soil_data_from_values, read_soil_data, soilparams
export read_management, validate_management, crop_inputs
export read_co2_series, climate_blocks, read_climate_block, climate_forcing

end
