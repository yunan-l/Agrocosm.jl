#!/usr/bin/env julia
include(joinpath(@__DIR__, "common.jl"))
length(ARGS) == 3 || error("usage: run_historical.jl STUDY.toml CFT_ID rainfed|irrigated")
config_path = abspath(ARGS[1])
config = TOML.parsefile(config_path)
root = check_study_config(config)
cft_id, water = parse(Int, ARGS[2]), ARGS[3]
patch = check_patch(config, cft_id, water)
haskey(ENV, "SLURM_JOB_ID") || error("global historical simulation is server/Slurm-only")

using MPI
include(joinpath(@__DIR__, "historical_helpers.jl"))
check_loaded_attribution_model()
MPI.Init()
try
    comm = MPI.COMM_WORLD
    rank, ranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    directory = joinpath(root, "historical_1901_2019", "batches", patch)
    rank == 0 && claim_stage_directory(config, "historical_1901_2019", patch)
    MPI.Barrier(comm)
    rank_directory = require_inside(joinpath(directory, "rank_$(lpad(rank, 4, '0'))"), root)
    mkdir(rank_directory)
    context = default_historical_context(config, config_path, cft_id, water; rank, ranks)
    simulation = context.simulation
    variables = production_output_variables()
    chunks = OutputChunk[]
    writer = chunk -> push!(chunks, chunk)
    stream = OutputStream(variables; frequency = :annual, writer,
        cell_ids = context.selection.cell_ids)
    for index in eachindex(context.years)
        run_output_block!(simulation, context.forcings[index], context.management_blocks[index], stream)
    end
    finish_output_stream!(stream, simulation.simulated_days)
    simulation.simulated_days == 365 * 119 || error("historical simulation did not finish")
    length(chunks) == 119 || error("historical output must contain all 119 years")
    output_path = joinpath(rank_directory, "annual_1901_2019.nc")
    # Keep each rank compact: do not allocate a global (lon, lat, 119) cube
    # independently on every MPI rank.
    write_compact_history(output_path, context.grid, context.selection, chunks, context.years)
    record = completion_record(config, config_path;
        stage = "historical_1901_2019", cft_id, water_system = water,
        partition_rank = rank, partition_count = ranks,
        cell_ids = context.selection.cell_ids, output_path, output_sha256 = file_sha(output_path),
        simulation_start_year = 1901, simulation_end_year = 2019,
        allocation_path = context.allocation_file,
        allocation_sha256 = file_sha(context.allocation_file))
    write_study_toml(joinpath(rank_directory, "study_complete.toml"), record, root)
    MPI.Barrier(comm)
    if rank == 0
        manifests = [joinpath(directory, "rank_$(lpad(index, 4, '0'))", "study_complete.toml")
            for index in 0:(ranks - 1)]
        records = TOML.parsefile.(manifests)
        ids = reduce(vcat, [item["cell_ids"] for item in records])
        allunique(ids) || error("historical MPI partitions overlap")
        allocation = read_soil_pool_allocation(context.allocation_file)
        Set(ids) == Set(allocation.selection.cell_ids) || error("historical cells do not cover the warm-up")
        all(item -> item["simulation_config_sha256"] == simulation_config_sha(config), records) ||
            error("historical rank configuration mismatch")
        write_study_toml(joinpath(directory, "study_complete.toml"),
            completion_record(config, config_path; stage = "historical_1901_2019",
                cft_id, water_system = water, partition_count = ranks,
                rank_manifests = manifests, cell_count = length(ids)), root)
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
