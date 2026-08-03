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

function run_global_cfts(config_path; backend_override = :cpu)
    config = TOML.parsefile(config_path)
    haskey(config["paths"], "pool_allocation") && error(
        "the CFT workflow derives one CFT- and water-system-specific allocation per batch; " *
        "do not set paths.pool_allocation in the shared configuration",
    )
    systems = requested_crop_systems(config)
    output_directory = abspath(config["paths"]["output_directory"])
    haskey(config, "free_warmup") && haskey(config, "allocation_validation") && error(
        "use [free_warmup], not both [free_warmup] and legacy [allocation_validation]",
    )
    free_warmup_options = get(config, "free_warmup",
        get(config, "allocation_validation", nothing))
    simulation_years = Int(config["run"]["simulation_start_year"]):Int(
        config["run"]["simulation_end_year"],
    )
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
        production_output = joinpath(
            result.output_directory,
            "global_wheat_$(first(simulation_years))_$(last(simulation_years)).nc",
        )
        push!(batch_manifest, Dict(
            "name" => name,
            "cft_id" => cft_id,
            "irrigated" => irrigated,
            "calibration_output_directory" => calibration.output_directory,
            "pool_allocation" => allocation_path,
            "production_output_directory" => result.output_directory,
            "production_output" => production_output,
            "production_warmup_years" => result.warmup_years,
            "production_warmup_converged" => result.warmup_converged,
        ))
    end
    manifest_path = joinpath(output_directory, "cft_batch_manifest.toml")
    write_report(manifest_path, Dict("batches" => batch_manifest))
    return manifest_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_cfts_cpu.jl CONFIG_TOML")
    println(run_global_cfts(abspath(ARGS[1]); backend_override = :cpu))
end
