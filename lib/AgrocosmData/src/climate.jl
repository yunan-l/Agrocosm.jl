const _CLIMATE_DATASETS = (:temp, :prec, :lwnet, :swdown)

"""Lazy, restartable iterator over bounded daily climate blocks."""
struct ClimateBlockReader{T <: AbstractFloat, TT, G, S, C}
    catalog::C
    grid::G
    selection::S
    time::Vector{TT}
    indices::UnitRange{Int}
    block_days::Int
    co2::CO2Series{T}
end

"""Read a two-column `year ppm` global annual CO₂ text file."""
function read_co2_series(path::AbstractString; T::Type{<:AbstractFloat} = Float32)
    years = Int32[]
    values = T[]
    for (line_number, line) in enumerate(eachline(path))
        content = strip(first(split(line, '#'; limit = 2)))
        isempty(content) && continue
        fields = split(content)
        length(fields) >= 2 || throw(ArgumentError("invalid CO₂ row $line_number in $path"))
        year = tryparse(Int32, fields[1])
        value = tryparse(T, fields[2])
        isnothing(year) && throw(ArgumentError("invalid CO₂ year on row $line_number in $path"))
        isnothing(value) && throw(ArgumentError("invalid CO₂ value on row $line_number in $path"))
        value > zero(T) || throw(ArgumentError("CO₂ values must be positive"))
        push!(years, year)
        push!(values, value)
    end
    isempty(years) && throw(ArgumentError("CO₂ file contains no data"))
    issorted(years) || throw(ArgumentError("CO₂ years must be sorted"))
    allunique(years) || throw(ArgumentError("CO₂ years must be unique"))
    provenance = DataProvenance(DATA_SCHEMA_VERSION, abspath(path), "co2", "ppm")
    return CO2Series(years, values, provenance)
end

function _climate_time(catalog::DatasetCatalog)
    reference = _time_coordinate(dataset(catalog, :temp), Colon())
    isempty(reference) && throw(ArgumentError("climate time coordinate cannot be empty"))
    for name in _CLIMATE_DATASETS[2:end]
        candidate = _time_coordinate(dataset(catalog, name), Colon())
        candidate == reference || throw(ArgumentError("$name time coordinate does not match temp"))
    end
    return reference
end

function _climate_indices(time, start_year, end_year)
    if isnothing(start_year) && isnothing(end_year)
        return 1:length(time)
    end
    first_year = isnothing(start_year) ? minimum(_calendar_year.(time)) : Int(start_year)
    last_year = isnothing(end_year) ? maximum(_calendar_year.(time)) : Int(end_year)
    first_year <= last_year || throw(ArgumentError("start_year must not exceed end_year"))
    indices = findall(value -> first_year <= _calendar_year(value) <= last_year, time)
    isempty(indices) && throw(ArgumentError("requested climate years are unavailable"))
    indices == collect(first(indices):last(indices)) ||
        throw(ArgumentError("requested climate rows must be contiguous"))
    return first(indices):last(indices)
end

"""Create a lazy iterator that reads at most `block_days` daily rows at once."""
function climate_blocks(
    catalog::DatasetCatalog,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    co2_path::AbstractString = dataset(catalog, :co2).path,
    start_year = nothing,
    end_year = nothing,
    block_days::Integer = 31,
    T::Type{<:AbstractFloat} = Float32,
)
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    time = _climate_time(catalog)
    indices = _climate_indices(time, start_year, end_year)
    co2 = read_co2_series(co2_path; T)
    requested_years = unique(Int32.(_calendar_year.(time[indices])))
    available_years = Set(co2.years)
    all(year -> year in available_years, requested_years) ||
        throw(ArgumentError("CO₂ data do not cover all requested climate years"))
    return ClimateBlockReader{T, eltype(time), typeof(grid), typeof(selection), typeof(catalog)}(
        catalog, grid, selection, time, indices, Int(block_days), co2,
    )
end

Base.length(reader::ClimateBlockReader) = cld(length(reader.indices), reader.block_days)

function _block_indices(reader::ClimateBlockReader, block_index::Integer)
    1 <= block_index <= length(reader) || throw(BoundsError(reader, block_index))
    first_index = first(reader.indices) + (block_index - 1) * reader.block_days
    last_index = min(last(reader.indices), first_index + reader.block_days - 1)
    return first_index:last_index
end

function _validate_climate(name::Symbol, values)
    any(ismissing, values) && throw(ArgumentError("$name contains missing values"))
    all(isfinite, values) || throw(ArgumentError("$name contains non-finite values"))
    name === :prec && !all(>=(0), values) && throw(ArgumentError("precipitation must be non-negative"))
    name === :swdown && !all(>=(0), values) && throw(ArgumentError("shortwave radiation must be non-negative"))
    name === :temp && !all(value -> -100 <= value <= 70, values) &&
        throw(ArgumentError("temperature must be supplied in degrees Celsius"))
    return values
end

function _daily_co2(series::CO2Series{T}, time) where {T}
    positions = Dict(year => index for (index, year) in pairs(series.years))
    return T[series.values[positions[Int32(_calendar_year(value))]] for value in time]
end

"""Read one numbered block from a `ClimateBlockReader`."""
function read_climate_block(reader::ClimateBlockReader{T}, block_index::Integer) where {T}
    indices = _block_indices(reader, block_index)
    fields = map(_CLIMATE_DATASETS) do name
        read_compact_variable(
            dataset(reader.catalog, name), reader.grid;
            selection = reader.selection,
            selectors = (time = indices,),
            order = (:time, :cell),
            T,
        )
    end
    temp, prec, lwnet, swdown = fields
    temp_values = Matrix{T}(_validate_climate(:temp, temp.values))
    prec_values = Matrix{T}(_validate_climate(:prec, prec.values))
    lwnet_values = Matrix{T}(_validate_climate(:lwnet, lwnet.values))
    swdown_values = Matrix{T}(_validate_climate(:swdown, swdown.values))
    block_time = collect(reader.time[indices])
    provenance = (
        temp = temp.provenance,
        prec = prec.provenance,
        lwnet = lwnet.provenance,
        swdown = swdown.provenance,
        co2 = reader.co2.provenance,
    )
    return ClimateBlock(
        block_time,
        temp_values,
        prec_values,
        swdown_values,
        lwnet_values,
        _daily_co2(reader.co2, block_time),
        reader.selection,
        provenance,
    )
end

function Base.iterate(reader::ClimateBlockReader, block_index::Int = 1)
    block_index > length(reader) && return nothing
    return read_climate_block(reader, block_index), block_index + 1
end

"""Return the model-facing forcing tuple without provenance metadata."""
climate_forcing(block::ClimateBlock) = (
    temp = block.temperature,
    prec = block.precipitation,
    sw = block.shortwave,
    lw = block.longwave,
    co2 = block.co2,
    co2_daily = true,
)
