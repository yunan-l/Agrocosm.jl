#!/usr/bin/env julia
include(joinpath(@__DIR__, "common.jl"))
length(ARGS) == 3 || error("usage: run_warmup.jl STUDY.toml CFT_ID rainfed|irrigated")
config_path = abspath(ARGS[1])
config = TOML.parsefile(config_path)
root = check_study_config(config)
cft_id, water = parse(Int, ARGS[2]), ARGS[3]
patch = check_patch(config, cft_id, water)
haskey(ENV, "SLURM_JOB_ID") || error("global warm-up is server/Slurm-only")

include(joinpath(@__DIR__, "..", "run_global_cfts_mpi.jl"))
check_loaded_attribution_model()
MPI.Init()
try
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    directory = joinpath(root, "warmup_600y", "batches", patch)
    generated = joinpath(directory, "warmup_config.toml")
    if rank == 0
        claim_stage_directory(config, "warmup_600y", patch)
        stage = deepcopy(config)
        stage["paths"]["output_directory"] = joinpath(root, "warmup_600y")
        # Retain rank-local convergence/selection audits for scientific QC.
        # Calibration-only ranks do not write full historical climate output.
        stage["mpi"] = Dict("cleanup_rank_outputs" => false)
        merge!(stage["run"], Dict(
            "calibration_only" => true, "resume_completed_batches" => false,
            "warmup_climate_start_year" => 1901, "warmup_climate_end_year" => 1930,
            "warmup_target_constrained" => true, "warmup_minimum_years" => 600,
            "warmup_maximum_years" => 600, "warmup_consecutive_years" => 3,
            "warmup_relative_tolerance" => 0.01, "warmup_pool_fraction_tolerance" => 0.01,
            "warmup_required_converged_fraction" => 1.0, "warmup_cache_climate" => false,
            "require_warmup_convergence" => false, "diagnostic_cells" => 10, "cell_limit" => 0,
        ))
        write_study_toml(generated, stage, root)
    end
    MPI.Barrier(comm)
    manifest = run_global_cfts_mpi(generated; cft_id, irrigated = water == "irrigated")
    if rank == 0
        product = only(TOML.parsefile(manifest)["batches"])
        allocation = require_inside(product["pool_allocation"], joinpath(root, "warmup_600y"))
        record = completion_record(config, config_path;
            stage = "warmup_600y", cft_id, water_system = water, warmup_years = 600,
            allocation_path = allocation, allocation_sha256 = file_sha(allocation),
            source_manifest = manifest)
        write_study_toml(joinpath(directory, "study_complete.toml"), record, root)
    end
    MPI.Barrier(comm)
catch exception
    showerror(stderr, exception, catch_backtrace())
    println(stderr)
    flush(stderr)
    MPI.Abort(MPI.COMM_WORLD, 1)
    rethrow()
finally
    MPI.Finalize()
end
