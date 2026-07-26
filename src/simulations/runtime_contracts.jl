abstract type AbstractExecutionArchitecture end

"""Host execution using ordinary Julia arrays."""
struct HostArchitecture <: AbstractExecutionArchitecture end

"""Array-backed accelerator execution, parameterized by its transfer function."""
struct AcceleratorArchitecture{F} <: AbstractExecutionArchitecture
    device::F
end

"""Compact active-cell domain and its stable external cell identifiers."""
struct ActiveLandDomain{I <: AbstractVector{Int}, C <: AbstractVector{Int32}}
    indices::I
    cell_ids::C

    function ActiveLandDomain(
        indices::AbstractVector{<:Integer},
        cell_ids::AbstractVector{<:Integer} = indices,
    )
        length(indices) == length(cell_ids) || throw(DimensionMismatch(
            "active indices and cell ids must have equal length",
        ))
        isempty(indices) && throw(ArgumentError("active land domain cannot be empty"))
        all(>(0), indices) || throw(ArgumentError("active indices must be positive"))
        allunique(indices) || throw(ArgumentError("active indices must be unique"))
        allunique(cell_ids) || throw(ArgumentError("active cell ids must be unique"))
        return new{Vector{Int}, Vector{Int32}}(Int.(indices), Int32.(cell_ids))
    end
end

"""Precision, architecture, and active-domain contract for one simulation."""
struct ExecutionContext{T <: AbstractFloat, A <: AbstractExecutionArchitecture, D}
    architecture::A
    domain::D
end

function ExecutionContext(
    ::Type{T}, device, indices::AbstractVector{<:Integer}; cell_ids = indices,
) where {T <: AbstractFloat}
    architecture = device === identity ? HostArchitecture() : AcceleratorArchitecture(device)
    domain = ActiveLandDomain(indices, cell_ids)
    return ExecutionContext{T, typeof(architecture), typeof(domain)}(architecture, domain)
end

float_type(::ExecutionContext{T}) where {T} = T
array_device(::HostArchitecture) = identity
array_device(architecture::AcceleratorArchitecture) = architecture.device
array_device(context::ExecutionContext) = array_device(context.architecture)
architecture_name(::HostArchitecture) = :cpu
architecture_name(::AcceleratorArchitecture) = :accelerator
architecture_name(context::ExecutionContext) = architecture_name(context.architecture)

function Base.show(io::IO, domain::ActiveLandDomain)
    print(io, "ActiveLandDomain(", length(domain.indices), " cells)")
end

function Base.show(io::IO, context::ExecutionContext)
    print(
        io,
        "ExecutionContext(", architecture_name(context), ", ",
        float_type(context), ", ", length(context.domain.indices), " cells)",
    )
end

"""Machine-readable metadata for a model or output variable."""
struct VariableSpec
    path::Tuple{Vararg{Symbol}}
    role::Symbol
    dimensions::Tuple{Vararg{Symbol}}
    units::String
    description::String
end

const _OUTPUT_VARIABLE_METADATA = Dict{Tuple{Symbol, Symbol}, NamedTuple}(
    (:crop, :gpp) => (units = "gC m-2 day-1", description = "Gross primary production"),
    (:crop, :npp) => (units = "gC m-2 day-1", description = "Net primary production"),
    (:crop, :lambda) => (units = "1", description = "Intercellular-to-ambient CO2 ratio"),
    (:crop, :potential_vcmax) => (units = "gC m-2 day-1", description = "Potential maximum carboxylation capacity"),
    (:crop, :vcmax) => (units = "gC m-2 day-1", description = "Realized maximum carboxylation capacity"),
    (:crop, :nitrogen_limitation) => (units = "1", description = "Nitrogen limitation factor"),
    (:crop, :respiration) => (units = "gC m-2 day-1", description = "Plant respiration"),
    (:crop, :biomass) => (units = "gC m-2", description = "Live crop carbon"),
    (:crop, :lai) => (units = "m2 m-2", description = "Leaf area index"),
    (:crop, :storage_carbon) => (units = "gC m-2", description = "Storage-organ carbon"),
    (:crop, :yield) => (units = "gC m-2 year-1", description = "Harvested storage-organ carbon"),
    (:crop, :fphu) => (units = "1", description = "Fraction of potential heat units"),
    (:crop, :water_deficit) => (units = "%", description = "Crop water deficit"),
    (:crop, :growing_mask) => (units = "1", description = "Active crop-stand mask"),
    (:calendar, :harvesting_mask) => (units = "1", description = "Harvest condition mask"),
    (:calendar, :harvesting_year) => (units = "year", description = "Simulation harvest year"),
    (:calendar, :harvest_date) => (units = "day_of_year", description = "Harvest day of year"),
    (:calendar, :sowing_event) => (units = "1", description = "Daily sowing event"),
    (:calendar, :harvest_event) => (units = "1", description = "Daily harvest event"),
)

function output_variable_spec(group::Symbol, field::Symbol)
    metadata = get(_OUTPUT_VARIABLE_METADATA, (group, field), nothing)
    isnothing(metadata) && throw(ArgumentError("missing metadata for output variable $group.$field"))
    frequency = field in (:yield, :harvest_date, :harvesting_year) ? :annual : :daily
    return VariableSpec(
        (:output, group, field), :output, (:time, :cell),
        metadata.units, metadata.description,
    ), frequency
end

function _state_array_specs!(specs, value, path::Tuple, role::Symbol, cells::Int)
    if value isa AbstractArray
        dimensions = ntuple(
            index -> index == ndims(value) && size(value, index) == cells ?
                :cell : Symbol(:axis_, index),
            ndims(value),
        )
        push!(specs, VariableSpec(path, role, dimensions, "", ""))
        return specs
    end
    value isa NamedTuple || isstructtype(typeof(value)) || return specs
    for name in propertynames(value)
        _state_array_specs!(specs, getproperty(value, name), (path..., name), role, cells)
    end
    return specs
end

"""Return the numerical state inventory grouped by lifecycle role."""
function state_schema(state::ModelState)
    specs = VariableSpec[]
    cells = length(state.inputs.weather.temp)
    for role in (:prognostic, :fluxes, :auxiliary, :inputs, :events, :workspace)
        _state_array_specs!(specs, getproperty(state, role), (role,), role, cells)
    end
    paths = getfield.(specs, :path)
    allunique(paths) || throw(ArgumentError("model state contains duplicate variable paths"))
    return specs
end

"""Validate that every runtime state array ends in the active-cell dimension."""
function validate_state_schema(state::ModelState, cells::Integer)
    cells > 0 || throw(ArgumentError("cells must be positive"))
    for spec in state_schema(state)
        value = foldl(getproperty, spec.path; init = state)
        :cell in spec.dimensions || continue
        size(value, ndims(value)) == cells || throw(DimensionMismatch(
            "$(join(spec.path, '.')) ends in size $(size(value, ndims(value))); expected $cells active cells",
        ))
    end
    return state
end
