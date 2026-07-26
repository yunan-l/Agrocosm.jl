using Dates
using NCDatasets
using TOML

const PFT_DIMENSIONS = Set(("pft", "cft", "crop", "band"))
const TIME_DIMENSIONS = Set(("time", "year"))
const ANNUAL_MANAGEMENT = Set((
    "landuse",
    "fertilizer",
    "manure",
    "residue_fraction",
    "sowing_date",
    "phu",
))

dimension_kind(name) = lowercase(String(name)) in PFT_DIMENSIONS ? :pft :
    lowercase(String(name)) in TIME_DIMENSIONS ? :time : :other

attributes(value) = Dict(String(key) => item for (key, item) in pairs(value.attrib)
    if String(key) != "_NCProperties")

function copy_attributes!(destination, source)
    for (key, value) in attributes(source)
        destination.attrib[key] = value
    end
    return destination
end

calendar_year(value::Real) = round(Int, value)
calendar_year(value) = Dates.year(value)

function validate_365_day_years(time, years, label)
    year_values = calendar_year.(time)
    for year in years
        indices = findall(==(year), year_values)
        length(indices) == 365 || throw(ArgumentError(
            "$label year $year contains $(length(indices)) rows; expected exactly 365",
        ))
        if !isempty(indices) && !(time[first(indices)] isa Real)
            any(index -> (Dates.month(time[index]), Dates.day(time[index])) == (2, 29), indices) &&
                throw(ArgumentError("$label year $year contains 29 February"))
            (Dates.month(time[first(indices)]), Dates.day(time[first(indices)])) == (1, 1) ||
                throw(ArgumentError("$label year $year does not start on 1 January"))
            (Dates.month(time[last(indices)]), Dates.day(time[last(indices)])) == (12, 31) ||
                throw(ArgumentError("$label year $year does not end on 31 December"))
        end
    end
    return nothing
end

function indices_for_years(
    dataset,
    variable_name::AbstractString,
    years;
    require_365_days::Bool = true,
)
    dimensions = String.(dimnames(dataset[variable_name]))
    time_position = findfirst(name -> dimension_kind(name) === :time, dimensions)
    isnothing(time_position) && throw(ArgumentError("$variable_name has no time dimension"))
    time_name = dimensions[time_position]
    time = collect(dataset[time_name][:])
    indices = findall(value -> calendar_year(value) in years, time)
    found = unique(calendar_year.(time[indices]))
    found == collect(years) || throw(ArgumentError(
        "$variable_name contains years $found, expected $(collect(years))",
    ))
    require_365_days && validate_365_day_years(time, years, variable_name)
    return time_name, indices
end

function selected_time_values(path, variable_name, years)
    return NCDataset(path, "r") do dataset
        time_name, indices = indices_for_years(dataset, variable_name, years)
        collect(dataset[time_name][indices])
    end
end

function define_subset_variable!(destination, source, name, selections)
    source_variable = variable(source, name)
    dimensions = String.(dimnames(source_variable))
    output = defVar(
        destination, name, eltype(source_variable), Tuple(dimensions);
        attrib = attributes(source_variable),
    )
    source_indices = ntuple(position -> selections[dimensions[position]], length(dimensions))
    output_indices = ntuple(position -> 1:length(source_indices[position]), length(dimensions))
    output[output_indices...] = source_variable[source_indices...]
    return output
end

function copy_main_variable!(destination, source, name, selections; chunk_length::Integer)
    source_variable = variable(source, name)
    dimensions = String.(dimnames(source_variable))
    output = defVar(
        destination, name, eltype(source_variable), Tuple(dimensions);
        attrib = attributes(source_variable), deflatelevel = 1, shuffle = true,
    )
    time_position = findfirst(dimension -> dimension_kind(dimension) === :time, dimensions)
    if isnothing(time_position)
        source_indices = ntuple(position -> selections[dimensions[position]], length(dimensions))
        output_indices = ntuple(position -> 1:length(source_indices[position]), length(dimensions))
        output[output_indices...] = source_variable[source_indices...]
        return output
    end

    selected_time = selections[dimensions[time_position]]
    for first_output in 1:chunk_length:length(selected_time)
        last_output = min(length(selected_time), first_output + chunk_length - 1)
        output_time = first_output:last_output
        source_indices = ntuple(length(dimensions)) do position
            position == time_position ? selected_time[output_time] : selections[dimensions[position]]
        end
        output_indices = ntuple(length(dimensions)) do position
            position == time_position ? output_time : 1:length(source_indices[position])
        end
        output[output_indices...] = source_variable[source_indices...]
    end
    return output
end

function subset_netcdf(
    input_path::AbstractString,
    output_path::AbstractString,
    variable_name::AbstractString;
    pft_index::Union{Nothing, Integer} = nothing,
    years::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    require_365_days::Bool = true,
    chunk_length::Integer = 31,
)
    chunk_length > 0 || throw(ArgumentError("chunk_length must be positive"))
    mkpath(dirname(output_path))
    NCDataset(input_path, "r") do source
        haskey(source, variable_name) || throw(ArgumentError(
            "variable '$variable_name' not found in $input_path",
        ))
        dimensions = String.(dimnames(source[variable_name]))
        selections = Dict{String, Vector{Int}}()
        for (position, dimension) in pairs(dimensions)
            selections[dimension] = collect(1:size(source[variable_name], position))
        end
        if !isnothing(pft_index)
            pft_dimensions = filter(dimension -> dimension_kind(dimension) === :pft, dimensions)
            length(pft_dimensions) == 1 || throw(ArgumentError(
                "$variable_name must contain exactly one pft/cft/crop/band dimension",
            ))
            pft_dimension = only(pft_dimensions)
            1 <= pft_index <= length(selections[pft_dimension]) || throw(BoundsError())
            selections[pft_dimension] = [Int(pft_index)]
        end
        if !isnothing(years)
            time_dimension, time_indices = indices_for_years(
                source, variable_name, years; require_365_days,
            )
            selections[time_dimension] = time_indices
        end

        isfile(output_path) && throw(ArgumentError("output already exists: $output_path"))
        NCDataset(output_path, "c") do destination
            copy_attributes!(destination, source)
            destination.attrib["agrocosm_subset_source"] = abspath(input_path)
            destination.attrib["agrocosm_subset_variable"] = variable_name
            !isnothing(pft_index) &&
                (destination.attrib["agrocosm_rainfed_pft_index"] = Int32(pft_index))
            !isnothing(years) &&
                (destination.attrib["agrocosm_selected_years"] = join(years, ","))

            for dimension in dimensions
                defDim(destination, dimension, length(selections[dimension]))
            end
            for dimension in dimensions
                haskey(source, dimension) || continue
                coordinate_dimensions = String.(dimnames(source[dimension]))
                all(haskey(selections, name) for name in coordinate_dimensions) || continue
                define_subset_variable!(destination, source, dimension, selections)
            end
            copy_main_variable!(
                destination, source, variable_name, selections; chunk_length,
            )
        end
    end
    return output_path
end

function subset_co2(input_path, output_path, years)
    isfile(output_path) && throw(ArgumentError("output already exists: $output_path"))
    mkpath(dirname(output_path))
    selected = String[]
    for line in eachline(input_path)
        content = strip(first(split(line, '#'; limit = 2)))
        isempty(content) && continue
        fields = split(content)
        length(fields) >= 2 || continue
        year = tryparse(Int, fields[1])
        !isnothing(year) && year in years && push!(selected, line)
    end
    found = sort([parse(Int, first(split(line))) for line in selected])
    found == sort(collect(years)) || throw(ArgumentError(
        "CO₂ file contains years $found, expected $(sort(collect(years)))",
    ))
    open(output_path, "w") do output
        println(output, "# Extracted from $(abspath(input_path))")
        foreach(line -> println(output, line), selected)
    end
    return output_path
end

function resolve_path(config_path, value)
    path = String(value)
    return isabspath(path) ? path : normpath(joinpath(dirname(config_path), path))
end

function prepare_subset(config_path::AbstractString)
    config = TOML.parsefile(config_path)
    settings = config["subset"]
    output_directory = resolve_path(config_path, settings["output_directory"])
    pft_index = Int(get(settings, "rainfed_pft_index", 1))
    management_year = Int(get(settings, "management_year", 2015))
    climate_years = Int.(get(settings, "climate_years", [2015, 2016]))
    chunk_length = Int(get(settings, "chunk_length", 31))

    for name in sort!(collect(keys(config["management"])))
        spec = config["management"][name]
        subset_netcdf(
            resolve_path(config_path, spec["input"]),
            joinpath(output_directory, spec["output"]),
            spec["variable"];
            pft_index,
            years = name in ANNUAL_MANAGEMENT ? [management_year] : nothing,
            require_365_days = false,
            chunk_length,
        )
    end

    climate_names = sort!(collect(keys(config["climate"])))
    reference = config["climate"][first(climate_names)]
    reference_path = resolve_path(config_path, reference["input"])
    selected_years = climate_years
    reference_time = selected_time_values(
        reference_path, reference["variable"], selected_years,
    )
    for name in climate_names
        spec = config["climate"][name]
        input_path = resolve_path(config_path, spec["input"])
        selected_time_values(input_path, spec["variable"], selected_years) == reference_time ||
            throw(ArgumentError("climate time coordinate for $name does not match the reference"))
        subset_netcdf(
            input_path,
            joinpath(output_directory, spec["output"]),
            spec["variable"];
            years = selected_years, chunk_length,
        )
    end
    if haskey(config, "co2")
        subset_co2(
            resolve_path(config_path, config["co2"]["input"]),
            joinpath(output_directory, config["co2"]["output"]),
            selected_years,
        )
    end
    println(
        "Prepared rainfed PFT $pft_index, management year $management_year, " *
        "and climate years $(join(selected_years, ", ")) in $output_directory",
    )
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: julia --project=lib/AgrocosmData $(basename(@__FILE__)) CONFIG.toml")
    prepare_subset(abspath(only(ARGS)))
end
