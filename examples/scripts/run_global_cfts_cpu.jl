using NCDatasets
using TOML

include(joinpath(@__DIR__, "run_global_wheat_cpu.jl"))

function requested_crop_systems(config)
    cfts = get(config, "cfts", Dict{String, Any}())
    requested = get(cfts, "cft_ids", [1])
    cft_ids = requested isa AbstractString && lowercase(requested) == "all" ?
        collect(1:length(CFTS)) : Int.(requested)
    isempty(cft_ids) && error("cfts.cft_ids must not be empty")
    all(id -> 1 <= id <= length(CFTS), cft_ids) ||
        error("cfts.cft_ids must be in 1:$(length(CFTS))")
    allunique(cft_ids) || error("cfts.cft_ids must be unique")
    systems = Symbol.(lowercase.(String.(get(cfts, "water_systems", ["rainfed"]))))
    all(system -> system in (:rainfed, :irrigated), systems) || error(
        "cfts.water_systems must contain rainfed and/or irrigated",
    )
    allunique(systems) || error("cfts.water_systems must be unique")
    return [(cft_id, system === :irrigated) for cft_id in cft_ids for system in systems]
end

batch_name(cft_id, irrigated) = "cft_$(lpad(cft_id, 2, '0'))_" *
    (irrigated ? "irrigated" : "rainfed")

const _FREE_WARMUP_RUN_KEYS = (
    "warmup_minimum_years",
    "warmup_maximum_years",
    "warmup_consecutive_years",
    "warmup_relative_tolerance",
    "warmup_pool_fraction_tolerance",
    "warmup_required_converged_fraction",
    "warmup_cache_climate",
    "require_warmup_convergence",
)

const _CALIBRATION_RUN_OPTIONS = Dict{String, Any}(
    "warmup_minimum_years" => 600,
    "warmup_maximum_years" => 600,
    "warmup_target_constrained" => true,
    # Calibration always writes its allocation product. The free warm-up stage
    # is responsible for deciding whether the resulting state is production-ready.
    "require_warmup_convergence" => false,
)

function write_batch_config(path, config; output_directory, pool_allocation = nothing,
    production::Bool, target_constrained::Bool, warmup_options = nothing)
    batch = deepcopy(config)
    paths = batch["paths"]
    run = batch["run"]
    paths["output_directory"] = output_directory
    if isnothing(pool_allocation)
        pop!(paths, "pool_allocation", nothing)
    else
        paths["pool_allocation"] = pool_allocation
    end
    run["production"] = production
    run["warmup_target_constrained"] = target_constrained
    if !production
        merge!(run, _CALIBRATION_RUN_OPTIONS)
    elseif !isnothing(warmup_options)
        for key in _FREE_WARMUP_RUN_KEYS
            haskey(warmup_options, key) && (run[key] = warmup_options[key])
        end
    end
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, batch; sorted = true)
    end
    return path
end

function patch_landfrac(catalog, grid, systems, source_years)
    total = zeros(Float32, length(source_years), length(grid.cell_ids))
    fractions = Vector{Any}(undef, length(systems))
    for (index, (cft_id, irrigated)) in pairs(systems)
        landuse = read_management(
            catalog, :landuse, grid, cft_id;
            simulation_years = source_years, T = Float32, irrigated,
        )
        total .+= landuse.values
        fractions[index] = landuse
    end
    maximum(total) <= 1.0f0 + 1.0f-6 || error(
        "selected CFT patch land fractions exceed one in at least one grid cell",
    )
    return fractions
end

function write_cft_yield(path, batch_paths, systems, patch_fractions, grid)
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
            defVar(output, "cft_id", Int32, ("band",))[:] = Int32[first(system) for system in systems]
            defVar(output, "irrigated", Int8, ("band",))[:] = Int8[last(system) for system in systems]
            yield = defVar(output, "yield", Float32, ("longitude", "latitude", "band", "time"))
            yield.attrib["long_name"] = "unweighted crop yield by CFT and water system"
            yield.attrib["aggregation"] = "none"
            landfrac = defVar(output, "landfrac", Float32, ("longitude", "latitude", "band", "time"))
            landfrac.attrib["long_name"] = "crop land fraction by CFT and water system"
            landfrac_sum = defVar(output, "landfrac_sum", Float32, ("longitude", "latitude", "time"))
            landfrac_sum.attrib["long_name"] = "sum of simulated crop land fractions"
            total_fraction = zeros(Float32, length(longitude), length(latitude), length(time))
            for (band, batch_path) in pairs(batch_paths)
                NCDataset(batch_path, "r") do input
                    input["longitude"][:] == longitude || error("longitude differs across crop batches")
                    input["latitude"][:] == latitude || error("latitude differs across crop batches")
                    input["time"][:] == time || error("time differs across crop batches")
                    fraction_data = patch_fractions[band]
                    fraction_data.time == Int.(time) || error("landfrac years differ across crop batches")
                    batch_yield = input["crop_yield"][:, :, :]
                    for time_index in eachindex(time)
                        fraction = expand_to_grid(
                            vec(fraction_data.values[time_index, :]), grid; fill_value = 0.0f0,
                        )
                        valid = isfinite.(batch_yield[:, :, time_index])
                        fraction[.!valid] .= 0.0f0
                        landfrac[:, :, band, time_index] = fraction
                        total_fraction[:, :, time_index] .+= fraction
                    end
                    yield[:, :, band, :] = batch_yield
                end
            end
            maximum(total_fraction) <= 1.0f0 + 1.0f-6 || error(
                "valid CFT patch land fractions exceed one in the reconstructed output",
            )
            landfrac_sum[:, :, :] = total_fraction
        end
    end
    return path
end

function run_global_cfts(config_path; backend_override = :cpu)
    config = TOML.parsefile(config_path)
    haskey(config, "cfts") || error("multi-CFT runs require a [cfts] configuration section")
    haskey(config["paths"], "pool_allocation") && error(
        "multi-CFT workflow derives one CFT- and water-system-specific allocation per batch; " *
        "do not set paths.pool_allocation in the shared configuration",
    )
    systems = requested_crop_systems(config)
    output_directory = abspath(config["paths"]["output_directory"])
    haskey(config, "free_warmup") && haskey(config, "allocation_validation") && error(
        "use [free_warmup], not both [free_warmup] and legacy [allocation_validation]",
    )
    free_warmup_options = get(config, "free_warmup",
        get(config, "allocation_validation", nothing))
    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    simulation_years = Int(config["run"]["simulation_start_year"]):Int(
        config["run"]["simulation_end_year"],
    )
    patch_fractions = patch_landfrac(catalog, grid, systems, management_source_years(config, simulation_years))
    batch_paths = String[]
    batch_manifest = Dict{String, Any}[]
    for (cft_id, irrigated) in systems
        name = batch_name(cft_id, irrigated)
        println("running $name")
        batch_directory = joinpath(output_directory, "batches", name)
        calibration_directory = joinpath(batch_directory, "calibration")
        calibration_config = write_batch_config(
            joinpath(batch_directory, "calibration.toml"), config;
            output_directory = calibration_directory,
            production = false,
            target_constrained = true,
        )
        calibration = run_global_wheat(
            calibration_config;
            backend_override,
            cft_id,
            irrigated,
        )
        allocation_path = joinpath(calibration.output_directory, "warmup_soil_pool_allocation.nc")
        isfile(allocation_path) || error("$name did not write a soil-pool allocation")
        production_directory = joinpath(batch_directory, "production")
        production_config = write_batch_config(
            joinpath(batch_directory, "production.toml"), config;
            output_directory = production_directory,
            pool_allocation = allocation_path,
            production = true,
            target_constrained = false,
            warmup_options = free_warmup_options,
        )
        result = run_global_wheat(production_config; backend_override, cft_id, irrigated)
        push!(batch_paths, joinpath(
            result.output_directory,
            "global_wheat_$(first(simulation_years))_$(last(simulation_years)).nc",
        ))
        push!(batch_manifest, Dict(
            "name" => name,
            "cft_id" => cft_id,
            "irrigated" => irrigated,
            "calibration_output_directory" => calibration.output_directory,
            "pool_allocation" => allocation_path,
            "production_output_directory" => result.output_directory,
            "production_warmup_years" => result.warmup_years,
            "production_warmup_converged" => result.warmup_converged,
        ))
    end
    output_path = joinpath(
        output_directory,
        "global_cft_yield_$(first(simulation_years))_$(last(simulation_years)).nc",
    )
    write_report(joinpath(output_directory, "cft_batch_manifest.toml"), Dict("batches" => batch_manifest))
    return write_cft_yield(output_path, batch_paths, systems, patch_fractions, grid)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_cfts_cpu.jl CONFIG_TOML")
    println(run_global_cfts(abspath(ARGS[1]); backend_override = :cpu))
end
