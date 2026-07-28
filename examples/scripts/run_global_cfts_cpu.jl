using NCDatasets
using TOML

include(joinpath(@__DIR__, "run_global_wheat_cpu.jl"))

function requested_crop_systems(config)
    cfts = get(config, "cfts", Dict{String, Any}())
    pft_ids = Int.(get(cfts, "pft_ids", collect(1:12)))
    all(1 .<= pft_ids .<= length(CROP_PFTS)) || error("cfts.pft_ids must be in 1:12")
    allunique(pft_ids) || error("cfts.pft_ids must be unique")
    systems = Symbol.(lowercase.(String.(get(cfts, "water_systems", ["rainfed", "irrigated"]))))
    all(system -> system in (:rainfed, :irrigated), systems) || error(
        "cfts.water_systems must contain rainfed and/or irrigated",
    )
    allunique(systems) || error("cfts.water_systems must be unique")
    return [(pft_id, system === :irrigated) for pft_id in pft_ids for system in systems]
end

batch_name(pft_id, irrigated) = "cft_$(lpad(pft_id, 2, '0'))_" *
    (irrigated ? "irrigated" : "rainfed")

function validate_patch_landfrac(catalog, grid, systems, source_years)
    total = zeros(Float32, length(source_years), length(grid.cell_ids))
    for (pft_id, irrigated) in systems
        landuse = read_management(
            catalog, :landuse, grid, pft_id;
            simulation_years = source_years, T = Float32, irrigated,
        )
        total .+= landuse.values
    end
    maximum(total) <= 1.0f0 + 1.0f-6 || error(
        "selected CFT patch land fractions exceed one in at least one grid cell",
    )
    return total
end

function write_cft_yield(path, batch_paths, systems)
    first_path = first(batch_paths)
    NCDataset(first_path, "r") do reference
        longitude = reference["longitude"][:]
        latitude = reference["latitude"][:]
        time = reference["time"][:]
        NCDataset(path, "c") do output
            defDim(output, "longitude", length(longitude))
            defDim(output, "latitude", length(latitude))
            defDim(output, "band", length(systems))
            defDim(output, "time", length(time))
            defVar(output, "longitude", eltype(longitude), ("longitude",))[:] = longitude
            defVar(output, "latitude", eltype(latitude), ("latitude",))[:] = latitude
            defVar(output, "time", eltype(time), ("time",))[:] = time
            defVar(output, "band", Int32, ("band",))[:] = Int32.(1:length(systems))
            defVar(output, "pft_id", Int32, ("band",))[:] = Int32[first(system) for system in systems]
            defVar(output, "irrigated", Int8, ("band",))[:] = Int8[last(system) for system in systems]
            yield = defVar(output, "yield", Float32, ("longitude", "latitude", "band", "time"))
            yield.attrib["long_name"] = "unweighted crop yield by CFT and water system"
            yield.attrib["aggregation"] = "none"
            for (band, batch_path) in pairs(batch_paths)
                NCDataset(batch_path, "r") do input
                    input["longitude"][:] == longitude || error("longitude differs across crop batches")
                    input["latitude"][:] == latitude || error("latitude differs across crop batches")
                    input["time"][:] == time || error("time differs across crop batches")
                    yield[:, :, band, :] = input["crop_yield"][:, :, :]
                end
            end
        end
    end
    return path
end

function run_global_cfts(config_path; backend_override = :cpu)
    config = TOML.parsefile(config_path)
    systems = requested_crop_systems(config)
    output_directory = abspath(config["paths"]["output_directory"])
    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    simulation_years = Int(config["run"]["simulation_start_year"]):Int(
        config["run"]["simulation_end_year"],
    )
    validate_patch_landfrac(
        catalog, grid, systems, management_source_years(config, simulation_years),
    )
    batch_paths = String[]
    for (pft_id, irrigated) in systems
        name = batch_name(pft_id, irrigated)
        println("running $name")
        result = run_global_wheat(
            config_path;
            backend_override,
            pft_id,
            irrigated,
            output_subdirectory = joinpath("batches", name),
        )
        push!(batch_paths, joinpath(
            result.output_directory,
            "global_wheat_$(config["run"]["simulation_start_year"])_$(config["run"]["simulation_end_year"]).nc",
        ))
    end
    years = simulation_years
    output_path = joinpath(output_directory, "global_cft_yield_$(first(years))_$(last(years)).nc")
    return write_cft_yield(output_path, batch_paths, systems)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_cfts_cpu.jl CONFIG_TOML")
    println(run_global_cfts(abspath(ARGS[1])))
end
