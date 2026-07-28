const DATA_SCHEMA_VERSION = v"0.1.0"

"""Trace one prepared variable back to its source file and units."""
struct DataProvenance
    schema_version::VersionNumber
    source::String
    variable::String
    units::String
end

"""One crop's file-band positions for rainfed and irrigated management."""
struct ManagementBands
    rainfed::Vector{Int32}
    irrigated::Vector{Int32}

    function ManagementBands(rainfed::AbstractVector{<:Integer}, irrigated::AbstractVector{<:Integer})
        length(rainfed) == length(irrigated) ||
            throw(ArgumentError("rainfed and irrigated band maps must have equal length"))
        isempty(rainfed) && throw(ArgumentError("management band maps cannot be empty"))
        rainfed32 = Int32.(rainfed)
        irrigated32 = Int32.(irrigated)
        all(>(0), rainfed32) || throw(ArgumentError("management bands are one-based and must be positive"))
        all(>(0), irrigated32) || throw(ArgumentError("management bands are one-based and must be positive"))
        allunique(rainfed32) || throw(ArgumentError("rainfed management bands must be unique"))
        allunique(irrigated32) || throw(ArgumentError("irrigated management bands must be unique"))
        return new(rainfed32, irrigated32)
    end
end

"""Location and expected metadata for one external dataset variable."""
struct DatasetSpec
    path::String
    variable::String
    units::String
    pft_ids::Vector{Int32}
    management_bands::Union{Nothing, ManagementBands}
end

function DatasetSpec(
    path::AbstractString,
    variable::AbstractString;
    units::AbstractString = "",
    pft_ids::AbstractVector{<:Integer} = Int[],
    rainfed_bands::AbstractVector{<:Integer} = Int[],
    irrigated_bands::AbstractVector{<:Integer} = Int[],
)
    ids = Int32.(pft_ids)
    allunique(ids) || throw(ArgumentError("dataset PFT ids must be unique"))
    has_bands = !isempty(rainfed_bands) || !isempty(irrigated_bands)
    has_bands && !isempty(ids) && throw(ArgumentError("use either pft_ids or explicit management bands, not both"))
    bands = has_bands ? ManagementBands(rainfed_bands, irrigated_bands) : nothing
    return DatasetSpec(String(path), String(variable), String(units), ids, bands)
end

"""Stable mapping between external PFT ids, names, and array positions."""
struct PFTRegistry
    ids::Vector{Int32}
    names::Vector{String}

    function PFTRegistry(ids::AbstractVector{<:Integer}, names::AbstractVector{<:AbstractString})
        length(ids) == length(names) || throw(ArgumentError("PFT ids and names must have equal length"))
        isempty(ids) && throw(ArgumentError("PFT registry cannot be empty"))
        ids32 = Int32.(ids)
        names_string = String.(names)
        allunique(ids32) || throw(ArgumentError("PFT ids must be unique"))
        allunique(names_string) || throw(ArgumentError("PFT names must be unique"))
        return new(ids32, names_string)
    end
end

"""Configured data sources and the explicit PFT ordering used by their files."""
struct DatasetCatalog
    datasets::Dict{Symbol, DatasetSpec}
    pfts::PFTRegistry
end

"""Canonical compact representation of the configured longitude/latitude grid."""
struct GridIndex{T <: AbstractFloat}
    longitude::Vector{T}
    latitude::Vector{T}
    cellid::Matrix{Int32}
    cell_ids::Vector{Int32}
    longitude_indices::Vector{Int32}
    latitude_indices::Vector{Int32}
end

"""Ordered subset of canonical compact cells."""
struct CellSelection
    compact_indices::Vector{Int}
    cell_ids::Vector{Int32}

    function CellSelection(compact_indices::AbstractVector{<:Integer}, cell_ids::AbstractVector{<:Integer})
        length(compact_indices) == length(cell_ids) ||
            throw(ArgumentError("compact indices and cell ids must have equal length"))
        return new(Int.(compact_indices), Int32.(cell_ids))
    end
end

"""One or more crop-system patches, allowing the same grid cell in multiple patches."""
struct PatchDomain{T <: AbstractFloat}
    patch_ids::Vector{Int32}
    compact_indices::Vector{Int}
    cell_ids::Vector{Int32}
    pft_ids::Vector{Int32}
    irrigated::BitVector
    landfrac::Vector{T}

    function PatchDomain(
        patch_ids::AbstractVector{<:Integer},
        compact_indices::AbstractVector{<:Integer},
        cell_ids::AbstractVector{<:Integer},
        pft_ids::AbstractVector{<:Integer},
        irrigated::AbstractVector{Bool},
        landfrac::AbstractVector{T},
    ) where {T <: AbstractFloat}
        lengths = length.((patch_ids, compact_indices, cell_ids, pft_ids, irrigated, landfrac))
        all(==(first(lengths)), lengths) || throw(DimensionMismatch(
            "every PatchDomain field must have the same length",
        ))
        patch_ids32 = Int32.(patch_ids)
        compact = Int.(compact_indices)
        pft_ids32 = Int32.(pft_ids)
        allunique(patch_ids32) || throw(ArgumentError("patch ids must be unique"))
        all(>(0), compact) || throw(ArgumentError("compact indices must be positive"))
        all(>(0), pft_ids32) || throw(ArgumentError("PFT ids must be positive"))
        all(isfinite, landfrac) && all(>=(zero(T)), landfrac) || throw(ArgumentError(
            "patch land fractions must be finite and non-negative",
        ))
        return new{T}(
            patch_ids32,
            compact,
            Int32.(cell_ids),
            pft_ids32,
            BitVector(irrigated),
            collect(landfrac),
        )
    end
end

"""An array with named canonical dimensions ending in `:cell`."""
struct CompactVariable{T, N, A <: AbstractArray{T, N}}
    values::A
    dimensions::NTuple{N, Symbol}
    selection::CellSelection
    provenance::DataProvenance
end

"""One selected PFT's time-varying management field."""
struct TimeCellData{T, TT, A <: AbstractMatrix{T}, V <: AbstractVector{TT}}
    time::V
    values::A
    selection::CellSelection
    pft_id::Int32
    irrigated::Bool
    provenance::DataProvenance
end

"""Annual global atmospheric CO₂ series in ppm."""
struct CO2Series{T <: AbstractFloat}
    years::Vector{Int32}
    values::Vector{T}
    provenance::DataProvenance
end

"""One bounded daily climate block aligned to a fixed compact cell selection."""
struct ClimateBlock{T <: AbstractFloat, TT}
    time::Vector{TT}
    temperature::Matrix{T}
    precipitation::Matrix{T}
    shortwave::Matrix{T}
    longwave::Matrix{T}
    co2::Vector{T}
    selection::CellSelection
    provenance::NamedTuple
end

"""Lazy model-facing view of climate blocks."""
struct ClimateForcingReader{R} <: AbstractVector{NamedTuple}
    source::R
end

"""Sequential forcing view that reads the following block on a worker thread."""
mutable struct PrefetchedClimateForcingReader{R} <: AbstractVector{NamedTuple}
    source::R
    scheduled_index::Int
    task::Union{Nothing, Task}
    closed::Bool
end

"""Fixed allocation selection and annual activity for one PFT."""
struct CropMask{T, A <: AbstractMatrix{T}}
    selection::CellSelection
    fraction::A
    active::BitMatrix
end

"""Versioned soil-code lookup used by the current Agrocosm parameterization."""
struct SoilLookup{T <: AbstractFloat}
    version::VersionNumber
    sand::Vector{T}
    silt::Vector{T}
    clay::Vector{T}
    saturation::Vector{T}
    diffusivity_dry::Vector{T}
    diffusivity_15::Vector{T}
    layer_depth::Vector{T}
end

"""Compact static soil inputs compatible with Agrocosm initialization."""
struct SoilData{T <: AbstractFloat}
    selection::CellSelection
    soilcode::Vector{Int32}
    ph::Vector{T}
    saturation::Matrix{T}
    sand::Vector{T}
    silt::Vector{T}
    clay::Vector{T}
    diffusivity_dry::Vector{T}
    diffusivity_15::Vector{T}
    layer_depth::Vector{T}
    provenance::NamedTuple
end

"""Layer-resolved HWSD-derived SOC and total-N targets on compact Agrocosm cells."""
struct SoilCNTargets{T <: AbstractFloat}
    selection::CellSelection
    layer_bounds::Matrix{T}
    soil_organic_carbon::Matrix{T}
    total_nitrogen::Matrix{T}
    coverage::Matrix{T}
    uncertain::BitMatrix
    provenance::NamedTuple
end

"""Bounded-memory area accumulator used while streaming HWSD raster tiles."""
mutable struct SoilCNAggregator{T <: AbstractFloat}
    carbon_sum::Matrix{T}
    nitrogen_sum::Matrix{T}
    carbon_area::Matrix{T}
    nitrogen_area::Matrix{T}
    uncertain_area::Matrix{T}
    target_area::Vector{T}
end
