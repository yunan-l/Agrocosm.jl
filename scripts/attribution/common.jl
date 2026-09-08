using TOML
using SHA
using Dates

const ATTRIBUTION_REPOSITORY = normpath(joinpath(@__DIR__, "..", ".."))
const ATTRIBUTION_STAGES = ("warmup_600y", "historical_1901_2019", "attribution")

function resolved_study_path(path)
    path = abspath(path)
    islink(path) && !ispath(path) && error("broken symlink: $path")
    ispath(path) && return realpath(path)
    parent = dirname(path)
    parent == path && error("cannot resolve path $path")
    return joinpath(resolved_study_path(parent), basename(path))
end

function require_inside(path, root)
    resolved = resolved_study_path(path)
    relative = relpath(resolved, resolved_study_path(root))
    (relative == ".." || startswith(relative, "../") || isabspath(relative)) &&
        error("writable path escapes the isolated study: $path")
    return resolved
end

file_sha(path) = open(io -> bytes2hex(sha256(io)), path, "r")
patch_name(cft_id, water) = "cft_$(lpad(cft_id, 2, '0'))_$(water)"

function simulation_config_sha(config)
    # Case/reference selection and resource requests may be added after
    # screening the completed historical run; scientific inputs may not.
    io = IOBuffer()
    TOML.print(io, Dict(key => config[key] for key in
        ("study", "paths", "climate", "management", "run", "cfts")); sorted = true)
    return bytes2hex(sha256(take!(io)))
end

function check_study_config(config; check_repository = true)
    study = config["study"]
    root = resolved_study_path(study["results_root"])
    basename(root) == "default_parameter_attribution_v1" ||
        error("results must use the dedicated default_parameter_attribution_v1 directory")
    basename(dirname(root)) == "output" && basename(dirname(dirname(root))) == "paper1_global_cell" ||
        error("results must remain inside Paper 1's output directory")
    study["parameter_source"] == "version_defaults" || error("fitted parameters are not allowed")
    config["run"]["nitrogen_limit_vcmax"] === true || error("nitrogen limitation must remain enabled")
    config["run"]["crop_resp_fix"] === true || error("crop_resp_fix must remain enabled")
    config["run"]["backend"] == "cpu" || error("this server workflow is CPU/MPI")
    config["management"]["mode"] == "fixed" && config["management"]["fixed_year"] == 2015 ||
        error("this study retains prescribed 2015 management")
    get(config["management"], "sowing_mode", "prescribed_sdate") == "prescribed_sdate" ||
        error("dynamic sowing is outside the fixed-event attribution contract")
    config["management"]["fertilizer"] == "yes" &&
        config["management"]["manure"] === true &&
        config["management"]["with_tillage"] === true ||
        error("retain fertilizer, manure and tillage in all three stages")
    all(id -> id in (1, 2, 3, 9), config["cfts"]["cft_ids"]) || error("unsupported CFT")
    all(water -> water in ("rainfed", "irrigated"), config["cfts"]["water_systems"]) ||
        error("invalid water system")
    config["climate"]["wind_variable"] == "windspeed" || error("wind variable must be windspeed")
    window_days = get(get(config, "attribution", Dict()), "window_days", 14)
    window_days isa Integer && !(window_days isa Bool) && window_days > 0 ||
        error("attribution.window_days must be a positive integer")
    for key in ("output_directory", "allocation_root", "pool_allocation")
        haskey(config["paths"], key) && error("$key is derived inside this study, not inherited")
    end
    if check_repository
        realpath(study["repository"]) == realpath(ATTRIBUTION_REPOSITORY) ||
            error("runner is not inside the specified isolated repository")
        development_root = realpath(dirname(ATTRIBUTION_REPOSITORY))
        basename(development_root) == "attribution" ||
            error("use the independent attribution checkout, not a production checkout")
        commit = readchomp(`git -C $ATTRIBUTION_REPOSITORY rev-parse HEAD`)
        commit == study["expected_commit"] || error("pinned source revision mismatch")
        isempty(readchomp(`git -C $ATTRIBUTION_REPOSITORY status --porcelain`)) ||
            error("freeze/commit the isolated source before server execution")
        require_inside(first(DEPOT_PATH), joinpath(development_root, "environment"))
    end
    return root
end

function check_patch(config, cft_id, water)
    cft_id in config["cfts"]["cft_ids"] || error("CFT not selected in this study")
    water in config["cfts"]["water_systems"] || error("water system not selected in this study")
    return patch_name(cft_id, water)
end

function claim_stage_directory(config, stage, name)
    stage in ATTRIBUTION_STAGES || error("invalid attribution stage")
    occursin(r"^[a-z0-9_]+$", name) || error("invalid stage directory name")
    root = check_study_config(config)
    directory = require_inside(joinpath(root, stage, "batches", name), root)
    mkpath(dirname(directory))
    mkdir(directory) # Atomic claim; deliberately refuse an existing/partial run.
    return directory
end

function write_study_toml(path, values, root)
    path = require_inside(path, root)
    ispath(path) && error("refusing to overwrite $path")
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, values; sorted = true)
    end
    return path
end

function completion_record(config, config_path; kwargs...)
    return Dict{String, Any}(
        "status" => "complete", "study_config_sha256" => file_sha(config_path),
        "simulation_config_sha256" => simulation_config_sha(config),
        "repository_commit" => config["study"]["expected_commit"],
        "julia_version" => string(VERSION),
        "environment_project_sha256" => file_sha(Base.active_project()),
        "environment_manifest_sha256" => file_sha(joinpath(dirname(Base.active_project()), "Manifest.toml")),
        "parameter_source" => "version_defaults", "nitrogen_limit_vcmax" => true,
        "created_at" => string(now()), (string(key) => value for (key, value) in kwargs)...,
    )
end

function read_stage_completion(config, config_path, stage, cft_id, water)
    root = check_study_config(config)
    path = require_inside(joinpath(root, stage, "batches", patch_name(cft_id, water), "study_complete.toml"), root)
    record = TOML.parsefile(path)
    record["status"] == "complete" || error("upstream stage is incomplete")
    record["simulation_config_sha256"] == simulation_config_sha(config) || error("upstream scientific configuration mismatch")
    record["repository_commit"] == config["study"]["expected_commit"] || error("upstream source mismatch")
    record["parameter_source"] == "version_defaults" || error("upstream parameters are not defaults")
    record["julia_version"] == string(VERSION) || error("upstream Julia version mismatch")
    record["environment_project_sha256"] == file_sha(Base.active_project()) &&
        record["environment_manifest_sha256"] == file_sha(joinpath(dirname(Base.active_project()), "Manifest.toml")) ||
        error("keep the same resolved Julia environment across all three stages")
    return record
end

function check_loaded_attribution_model()
    require_inside(pathof(Agrocosm), ATTRIBUTION_REPOSITORY)
    return nothing
end
