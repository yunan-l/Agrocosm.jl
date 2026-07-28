function _variable_dimension_info(spec::DatasetSpec)
    return NCDataset(spec.path, "r") do ds
        haskey(ds, spec.variable) || throw(ArgumentError("variable '$(spec.variable)' not found in $(spec.path)"))
        variable = ds[spec.variable]
        source = Symbol.(dimnames(variable))
        canonical = _canonical_dimension.(source)
        return source, canonical, size(variable)
    end
end

function _time_coordinate(spec::DatasetSpec, selector)
    return NCDataset(spec.path, "r") do ds
        variable = ds[spec.variable]
        source_names = Symbol.(dimnames(variable))
        canonical_names = _canonical_dimension.(source_names)
        time_position = findfirst(==(:time), canonical_names)
        isnothing(time_position) && return [0]
        coordinate = _read_all(ds[String(source_names[time_position])])
        return collect(coordinate[_preserving_selector(selector)])
    end
end

function _time_metadata(spec::DatasetSpec)
    return NCDataset(spec.path, "r") do ds
        variable = ds[spec.variable]
        source_names = Symbol.(dimnames(variable))
        canonical_names = _canonical_dimension.(source_names)
        time_position = findfirst(==(:time), canonical_names)
        isnothing(time_position) && return (calendar = "", units = "")
        coordinate = ds[String(source_names[time_position])]
        calendar = haskey(coordinate.attrib, "calendar") ?
            lowercase(String(coordinate.attrib["calendar"])) : ""
        units = haskey(coordinate.attrib, "units") ? String(coordinate.attrib["units"]) : ""
        return (; calendar, units)
    end
end

function _calendar_year(value)
    value isa Integer && return Int(value)
    value isa AbstractFloat && return round(Int, value)
    return Dates.year(value)
end

function _indices_for_years(spec::DatasetSpec, years)
    isnothing(years) && return Colon()
    coordinate = _time_coordinate(spec, Colon())
    requested = Set(Int.(years))
    indices = findall(value -> _calendar_year(value) in requested, coordinate)
    isempty(indices) && throw(ArgumentError("none of the requested years occur in $(spec.path)"))
    found = Set(_calendar_year(coordinate[index]) for index in indices)
    found == requested || throw(ArgumentError("requested years $(sort!(collect(requested))) are not all available"))
    return indices
end

function _indices_for_simulation_years(spec::DatasetSpec, simulation_years)
    requested = Int.(collect(simulation_years))
    isempty(requested) && throw(ArgumentError("simulation_years must not be empty"))
    coordinate = _time_coordinate(spec, Colon())
    available = _calendar_year.(coordinate)
    length(unique(available)) == length(available) || throw(ArgumentError(
        "management time coordinate in $(spec.path) contains duplicate years",
    ))
    issorted(available) || throw(ArgumentError(
        "management years in $(spec.path) must be sorted",
    ))
    first_year, last_year = extrema(available)
    source_years = clamp.(requested, first_year, last_year)
    positions = Dict(year => index for (index, year) in pairs(available))
    all(haskey(positions, year) for year in source_years) || throw(ArgumentError(
        "management years in $(spec.path) are not continuous between $first_year and $last_year",
    ))
    indices = [positions[year] for year in source_years]
    unique_indices = unique(indices)
    lookup = Dict(index => position for (position, index) in pairs(unique_indices))
    return unique_indices, [lookup[index] for index in indices], requested
end

function _materialize_management(values, ::Type{T}; missing_value::Real = 0) where {T <: AbstractFloat}
    output = Matrix{T}(undef, size(values))
    for index in eachindex(values)
        value = values[index]
        output[index] = ismissing(value) ? T(missing_value) : T(value)
    end
    return output
end

function _active_values(values, active)
    isnothing(active) && return vec(values)
    size(active) == size(values) || throw(DimensionMismatch("active mask must match management values"))
    return values[active]
end

"""Validate physical ranges for one canonical management variable."""
function validate_management(name::Symbol, values::AbstractMatrix; active::Union{Nothing, AbstractMatrix{Bool}} = nothing)
    all(isfinite, values) || throw(ArgumentError("$name contains non-finite values"))
    checked = _active_values(values, active)
    if name in (:landuse, :landfrac, :residue_fraction, :residuefrac, :tillage, :with_tillage)
        all(value -> 0 <= value <= 1, checked) || throw(ArgumentError("$name must be in [0, 1]"))
    elseif name in (:sowing_date, :sdate)
        all(value -> 1 <= value <= 366, checked) || throw(ArgumentError("active sowing dates must be in 1:366"))
    elseif name in (:phu, :phusum)
        all(value -> value != 0, checked) || throw(ArgumentError(
            "active PHU values must be non-zero; negative values denote winter crops",
        ))
    elseif name in (:fertilizer, :manure)
        all(>=(0), checked) || throw(ArgumentError("$name must be non-negative"))
    end
    return values
end

"""
Read one management variable for a selected PFT as `time × cell`.

PFT-independent files are accepted and applied to the requested PFT. Missing
values are replaced with `missing_value`; callers should pass the land-use
activity mask to enforce active-year ranges for sowing date and PHU.
"""
function read_management(
    spec::DatasetSpec,
    name::Symbol,
    grid::GridIndex,
    registry::PFTRegistry,
    pft_id::Integer;
    selection::CellSelection = all_cells(grid),
    years = nothing,
    simulation_years = nothing,
    time_indices = nothing,
    active::Union{Nothing, AbstractMatrix{Bool}} = nothing,
    missing_value::Real = 0,
    T::Type{<:AbstractFloat} = Float32,
    irrigated::Bool = false,
)
    specified_time_selectors = count(!isnothing, (years, simulation_years, time_indices))
    specified_time_selectors <= 1 || throw(ArgumentError(
        "specify only one of years, simulation_years, or time_indices",
    ))
    source_names, canonical_names, variable_size = _variable_dimension_info(spec)
    time_position = findfirst(==(:time), canonical_names)
    pft_position = findfirst(==(:pft), canonical_names)
    file_pft_ids = isempty(spec.pft_ids) ? registry.ids : spec.pft_ids
    if !isnothing(pft_position) && spec.management_bands === nothing &&
            variable_size[pft_position] != length(file_pft_ids)
        throw(DimensionMismatch("PFT dimension has $(variable_size[pft_position]) entries but its configured mapping has $(length(file_pft_ids))"))
    end

    mapped_rows = nothing
    requested_years = nothing
    selected_time = if !isnothing(time_indices)
        time_indices
    elseif !isnothing(simulation_years) && !isnothing(time_position)
        indices, mapped_rows, requested_years = _indices_for_simulation_years(
            spec, simulation_years,
        )
        indices
    else
        _indices_for_years(spec, years)
    end
    selector_pairs = Pair{Symbol, Any}[]
    !isnothing(time_position) && push!(selector_pairs, :time => selected_time)
    if !isnothing(pft_position)
        crop_position = pft_index(registry, pft_id)
        file_band = if spec.management_bands === nothing
            position = findfirst(==(Int32(pft_id)), file_pft_ids)
            isnothing(position) && throw(ArgumentError("PFT id $pft_id is not present in $(spec.path)"))
            position
        else
            bands = irrigated ? spec.management_bands.irrigated : spec.management_bands.rainfed
            Int(bands[crop_position])
        end
        file_band <= variable_size[pft_position] ||
            throw(DimensionMismatch("configured band $file_band exceeds the file's $(variable_size[pft_position])-band PFT dimension"))
        push!(selector_pairs, :pft => file_band)
    end
    selectors = (; selector_pairs...)
    order = isnothing(time_position) ?
        (isnothing(pft_position) ? (:cell,) : (:pft, :cell)) :
        (isnothing(pft_position) ? (:time, :cell) : (:time, :pft, :cell))
    compact = read_compact_variable(spec, grid; selection, selectors, order)

    values = compact.values
    !isnothing(pft_position) && (values = dropdims(values; dims = findfirst(==(:pft), order)))
    isnothing(time_position) && (values = reshape(values, 1, :))
    if !isnothing(simulation_years)
        requested_years = Int.(collect(simulation_years))
        if isnothing(time_position)
            values = repeat(values; outer = (length(requested_years), 1))
        else
            values = values[mapped_rows, :]
        end
    end
    matrix = _materialize_management(values, T; missing_value)
    validate_management(name, matrix; active)
    time = if !isnothing(simulation_years)
        requested_years
    elseif isnothing(time_position)
        [0]
    else
        _time_coordinate(spec, selected_time)
    end
    return TimeCellData(time, matrix, selection, Int32(pft_id), irrigated, compact.provenance)
end

function read_management(
    catalog::DatasetCatalog,
    name::Symbol,
    grid::GridIndex,
    pft_id::Integer;
    kwargs...,
)
    return read_management(dataset(catalog, name), name, grid, catalog.pfts, pft_id; kwargs...)
end

function _aligned_management_values(name::Symbol, data::TimeCellData, reference::TimeCellData)
    data.selection.cell_ids == reference.selection.cell_ids ||
        throw(ArgumentError("$name uses a different cell selection"))
    data.pft_id == reference.pft_id || throw(ArgumentError("$name uses a different PFT"))
    data.irrigated == reference.irrigated ||
        throw(ArgumentError("$name uses a different rainfed/irrigated system"))
    data.time == reference.time || throw(ArgumentError("$name uses different management years"))
    return data.values
end

"""Build annual model-facing management rows aligned to simulation years."""
function management_schedule(;
    sowing_date::TimeCellData,
    phu::TimeCellData,
    manure::Union{Nothing, TimeCellData} = nothing,
    fertilizer::Union{Nothing, TimeCellData} = nothing,
    residue_fraction::TimeCellData,
    fertilizer_mode = :yes,
    manure_enabled::Bool = true,
)
    fertilizer_mode = fertilizer_mode isa Symbol ? fertilizer_mode :
        Symbol(lowercase(String(fertilizer_mode)))
    fertilizer_mode in (:no, :yes, :auto) ||
        throw(ArgumentError("fertilizer_mode must be :no, :yes, or :auto"))
    fertilizer_mode === :yes && fertilizer === nothing &&
        throw(ArgumentError("prescribed fertilization requires fertilizer data"))
    manure_enabled && manure === nothing &&
        throw(ArgumentError("prescribed manure requires manure data"))

    reference = sowing_date
    phu_values = _aligned_management_values(:phu, phu, reference)
    residue_values = _aligned_management_values(
        :residue_fraction, residue_fraction, reference,
    )
    zero_values = zeros(eltype(phu_values), size(phu_values))
    fertilizer_values = fertilizer_mode === :yes ?
        _aligned_management_values(:fertilizer, fertilizer, reference) : zero_values
    manure_values = manure_enabled ?
        _aligned_management_values(:manure, manure, reference) : zero_values
    return (
        years = Int32.(_calendar_year.(reference.time)),
        sdate = Int32.(round.(reference.values)),
        phu = phu_values,
        manure = manure_values,
        fertilizer = fertilizer_values,
        residuefrac = residue_values,
    )
end

function _single_management_row(name::Symbol, data::TimeCellData, selection, pft_id, irrigated)
    data.selection.cell_ids == selection.cell_ids || throw(ArgumentError("$name uses a different cell selection"))
    data.pft_id == pft_id || throw(ArgumentError("$name uses a different PFT"))
    data.irrigated == irrigated || throw(ArgumentError("$name uses a different rainfed/irrigated system"))
    size(data.values, 1) == 1 || throw(ArgumentError("$name must contain exactly one selected time"))
    return vec(data.values)
end

"""Build the current Agrocosm crop-input NamedTuple from single-time management data."""
function crop_inputs(;
    sowing_date::TimeCellData,
    phu::TimeCellData,
    manure::Union{Nothing, TimeCellData} = nothing,
    fertilizer::Union{Nothing, TimeCellData} = nothing,
    residue_fraction::TimeCellData,
    fertilizer_mode = :yes,
    manure_enabled::Bool = true,
)
    selection = sowing_date.selection
    pft_id = sowing_date.pft_id
    irrigated = sowing_date.irrigated
    sdate_values = _single_management_row(:sowing_date, sowing_date, selection, pft_id, irrigated)
    phu_values = _single_management_row(:phu, phu, selection, pft_id, irrigated)
    fertilizer_mode = fertilizer_mode isa Symbol ? fertilizer_mode : Symbol(lowercase(String(fertilizer_mode)))
    fertilizer_mode in (:no, :yes, :auto) ||
        throw(ArgumentError("fertilizer_mode must be :no, :yes, or :auto"))
    if fertilizer_mode === :yes && fertilizer === nothing
        throw(ArgumentError("prescribed fertilization requires fertilizer data"))
    end
    if manure_enabled && manure === nothing
        throw(ArgumentError("prescribed manure requires manure data"))
    end
    manure_values = !manure_enabled ? zeros(eltype(phu_values), length(phu_values)) :
        _single_management_row(:manure, manure, selection, pft_id, irrigated)
    fertilizer_values = if fertilizer_mode === :yes
        _single_management_row(:fertilizer, fertilizer, selection, pft_id, irrigated)
    else
        zeros(eltype(phu_values), length(phu_values))
    end
    residue_values = _single_management_row(:residue_fraction, residue_fraction, selection, pft_id, irrigated)
    return (
        sdate = Int32.(round.(sdate_values)),
        phu = phu_values,
        manure = manure_values,
        fertilizer = fertilizer_values,
        residuefrac = residue_values,
    )
end
