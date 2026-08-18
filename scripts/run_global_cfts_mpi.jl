using MPI

include(joinpath(@__DIR__, "run_global_cfts_cpu.jl"))

function run_global_cfts_mpi(
    config_path;
    cft_id::Union{Nothing, Integer} = nothing,
    irrigated::Union{Nothing, Bool} = nothing,
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    count = MPI.Comm_size(comm)
    config = TOML.parsefile(config_path)
    configured_output = abspath(config["paths"]["output_directory"])
    mpi_config = get(config, "mpi", Dict{String, Any}())
    output_base = abspath(get(mpi_config, "output_directory", configured_output))
    output_root = joinpath(output_base, "mpi_$(lpad(count, 4, '0'))_ranks")
    rank_directory = joinpath(output_root, "rank_$(lpad(rank, 4, '0'))")

    println("MPI rank $rank/$count: output=$rank_directory")
    rank_manifest = run_global_cfts(
        config_path;
        backend_override = :cpu,
        cft_id,
        irrigated,
        output_directory_override = rank_directory,
        partition_rank = rank,
        partition_count = count,
    )
    MPI.Barrier(comm)

    if rank == 0
        rank_directories = [joinpath(output_root, "rank_$(lpad(value, 4, '0'))")
                            for value in 0:(count - 1)]
        rank_manifests = [joinpath(path, "cft_batch_manifest.toml")
                          for path in rank_directories]
        all(isfile, rank_manifests) || error("one or more MPI rank manifests are missing")
        write_report(joinpath(output_root, "mpi_manifest.toml"), Dict(
            "schema_version" => "1",
            "repository_commit" => _repository_commit(),
            "config_path" => abspath(config_path),
            "config_fingerprint" => _file_fingerprint(config_path),
            "rank_count" => count,
            "rank_directories" => rank_directories,
            "rank_manifests" => rank_manifests,
            "created_at" => string(now()),
        ))
    end
    MPI.Barrier(comm)
    return rank_manifest
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 3) || error(
        "usage: run_global_cfts_mpi.jl CONFIG_TOML [CFT_ID rainfed|irrigated]",
    )
    cft_id = length(ARGS) == 3 ? parse(Int, ARGS[2]) : nothing
    irrigated = length(ARGS) == 3 ?
        Symbol(lowercase(ARGS[3])) === :irrigated : nothing
    length(ARGS) == 3 && Symbol(lowercase(ARGS[3])) in (:rainfed, :irrigated) ||
        length(ARGS) == 1 || error("water system must be rainfed or irrigated")

    MPI.Init()
    try
        println(run_global_cfts_mpi(abspath(ARGS[1]); cft_id, irrigated))
    catch exception
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        println(stderr, "MPI rank $rank failed: ", sprint(showerror, exception, catch_backtrace()))
        MPI.Abort(MPI.COMM_WORLD, 1)
    finally
        MPI.Finalize()
    end
end
