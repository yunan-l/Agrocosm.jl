#!/usr/bin/env julia

"""Run one production patch across MPI ranks and publish one merged output.

    mpiexecjl -n 32 julia --project=. scripts/run_global_production_mpi.jl \\
        CONFIG.toml CFT_ID rainfed|irrigated

Structure copied from the reference workflow's production wrapper
(`work/diff_task/paper1_global_cell/scripts/global_simulation/runtime/
run_global_historical_mpi.jl`): a resume check agreed by every rank, per-rank
directories under `.mpi_ranks/`, a barrier, then rank 0 merges the rank
manifests into one canonical-grid NetCDF and removes the shards only after both
final products exist.

This wrapper merges ONE product - the annual production NetCDF. The two-phase
wrapper also merges a `warmup_soil_pool_allocation.nc`, because its calibration
phase produces one; production reads that allocation as an input and produces
none.
"""

using MPI
using NCDatasets
using TOML
using Dates

include(joinpath(@__DIR__, "run_global_production_cpu.jl"))

"""Merge per-rank production outputs onto the canonical grid.

The invariant that matters: every selected cell is written by exactly one rank.
`partition_cell_ids` from each rank's manifest is what makes that checkable, and
an overlap is an error rather than a last-writer-wins blend - a blended cell
would look entirely plausible in the output.
"""
function merge_production_partitions(manifest_paths, output_path)
    manifests = TOML.parsefile.(manifest_paths)
    reference_manifest = first(manifests)
    for manifest in manifests
        for key in (
            "cft_id", "water_system", "processes_configuration", "rate_scale",
            "depletion_fraction",
            "simulation_start_year", "simulation_end_year",
            "crop_resp_fix", "nitrogen_limit_vcmax",
            "config_fingerprint", "allocation_path", "pool_allocation_cell_policy",
        )
            manifest[key] == reference_manifest[key] || error(
                "production MPI partitions differ in $key",
            )
        end
    end
    ranks = Int.(getindex.(manifests, "partition_rank"))
    count = Int(reference_manifest["partition_count"])
    count == length(manifests) || error("production MPI partition count mismatch")
    sort(ranks) == collect(0:(count - 1)) || error("production MPI ranks are incomplete")

    output_paths = String.(getindex.(manifests, "output_path"))
    all(isfile, output_paths) || error("one or more production MPI outputs are missing")
    reference = NCDataset(first(output_paths), "r")
    try
        longitude = Float32.(reference["longitude"][:])
        latitude = Float32.(reference["latitude"][:])
        time = Int32.(reference["time"][:])
        cellid = Int32.(reference["cellid"][:, :])
        coordinates = Set(("longitude", "latitude", "time", "cellid"))
        variable_names = sort!(filter(
            name -> !(name in coordinates), String.(collect(keys(reference))),
        ))
        expected_names = union(coordinates, Set(variable_names))
        positions = Dict{Int32, CartesianIndex{2}}(
            cellid[index] => index for index in CartesianIndices(cellid)
        )
        merged = Dict(name => fill(
            Float32(NaN), length(longitude), length(latitude), length(time),
        ) for name in variable_names)
        assigned = Set{Int32}()

        for (path, manifest) in zip(output_paths, manifests)
            selected = Int32.(manifest["partition_cell_ids"])
            indices = map(selected) do selected_cell
                selected_cell in assigned && error(
                    "production MPI partitions overlap at cell $selected_cell",
                )
                haskey(positions, selected_cell) || error(
                    "production MPI cell $selected_cell is absent from the canonical grid",
                )
                push!(assigned, selected_cell)
                return positions[selected_cell]
            end
            NCDataset(path, "r") do dataset
                Set(String.(collect(keys(dataset)))) == expected_names || error(
                    "production MPI output variables differ across ranks",
                )
                Float32.(dataset["longitude"][:]) == longitude ||
                    error("production MPI longitude coordinates differ")
                Float32.(dataset["latitude"][:]) == latitude ||
                    error("production MPI latitude coordinates differ")
                Int32.(dataset["time"][:]) == time ||
                    error("production MPI time coordinates differ")
                Int32.(dataset["cellid"][:, :]) == cellid ||
                    error("production MPI cell IDs differ")
                for name in variable_names
                    values = dataset[name][:, :, :]
                    for index in indices
                        merged[name][index[1], index[2], :] .= values[index[1], index[2], :]
                    end
                end
            end
        end

        mkpath(dirname(output_path))
        isfile(output_path) && rm(output_path; force = true)
        NCDataset(output_path, "c") do dataset
            defDim(dataset, "longitude", length(longitude))
            defDim(dataset, "latitude", length(latitude))
            defDim(dataset, "time", length(time))
            defVar(dataset, "longitude", Float32, ("longitude",))[:] = longitude
            defVar(dataset, "latitude", Float32, ("latitude",))[:] = latitude
            defVar(dataset, "time", Int32, ("time",))[:] = time
            defVar(dataset, "cellid", Int32, ("longitude", "latitude"))[:, :] = cellid
            for name in variable_names
                defVar(dataset, name, Float32, ("longitude", "latitude", "time"))[:, :, :] =
                    merged[name]
            end
            dataset.attrib["source"] = "Agrocosm production MPI partition merge"
            dataset.attrib["partition_count"] = count
        end
    finally
        close(reference)
    end
    return output_path
end

function run_global_production_mpi(args = ARGS)
    length(args) == 3 || error(
        "usage: run_global_production_mpi.jl CONFIG_TOML CFT_ID rainfed|irrigated",
    )
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    count = MPI.Comm_size(comm)
    config_path = abspath(args[1])
    cft_id = parse(Int, args[2])
    water = lowercase(args[3])
    water in ("rainfed", "irrigated") || error("invalid water system: $water")

    config = TOML.parsefile(config_path)
    output_root = abspath(config["paths"]["output_directory"])
    patch = "cft_$(lpad(cft_id, 2, '0'))_$(water)"
    final_directory = joinpath(output_root, "batches", patch, "production")
    final_manifest = joinpath(final_directory, "run_manifest.toml")

    # Resume by consensus: one rank finding a complete product is not enough,
    # because a partial re-run would then blend two versions.
    complete = false
    if isfile(final_manifest)
        manifest = TOML.parsefile(final_manifest)
        complete = get(manifest, "config_fingerprint", "") ==
                   _file_fingerprint(config_path) &&
                   Int(get(manifest, "partition_count", -1)) == count &&
                   isfile(String(get(manifest, "output_path", "")))
    end
    if MPI.Allreduce(complete ? Int32(1) : Int32(0), +, comm) == count
        rank == 0 && println("production MPI patch already complete: $final_manifest")
        return final_manifest
    end

    temporary_root = joinpath(
        output_root, ".mpi_ranks", patch, "mpi_$(lpad(count, 4, '0'))_ranks",
    )
    rank_directory = joinpath(temporary_root, "rank_$(lpad(rank, 4, '0'))")
    production_main([
        config_path, string(cft_id), water,
        string(rank), string(count), rank_directory,
    ])
    MPI.Barrier(comm)

    if rank == 0
        rank_directories = [
            joinpath(temporary_root, "rank_$(lpad(value, 4, '0'))")
            for value in 0:(count - 1)
        ]
        rank_manifests = [joinpath(directory, "run_manifest.toml")
                          for directory in rank_directories]
        all(isfile, rank_manifests) ||
            error("one or more production rank manifests are missing")
        first_rank = TOML.parsefile(first(rank_manifests))
        final_output = joinpath(final_directory, basename(String(first_rank["output_path"])))
        merge_production_partitions(rank_manifests, final_output)
        write_report(final_manifest, Dict(
            "schema_version" => "1",
            "entry" => "run_global_production_mpi.jl",
            "warmup" => "none: continued from the calibration allocation",
            "repository_commit" => _repository_commit(),
            "config_path" => config_path,
            "config_fingerprint" => _file_fingerprint(config_path),
            "cft_id" => cft_id,
            "water_system" => water,
            "processes_configuration" => first_rank["processes_configuration"],
            "rate_scale" => first_rank["rate_scale"],
            "depletion_fraction" => first_rank["depletion_fraction"],
            "management_mode" => first_rank["management_mode"],
            "management_fixed_year" => first_rank["management_fixed_year"],
            "crop_resp_fix" => first_rank["crop_resp_fix"],
            "nitrogen_limit_vcmax" => first_rank["nitrogen_limit_vcmax"],
            "allocation_root" => first_rank["allocation_root"],
            "allocation_path" => first_rank["allocation_path"],
            "pool_allocation_cell_policy" => first_rank["pool_allocation_cell_policy"],
            "candidate_cell_count" => first_rank["candidate_cell_count"],
            "allocation_excluded_cell_count" => first_rank["allocation_excluded_cell_count"],
            "simulation_start_year" => first_rank["simulation_start_year"],
            "simulation_end_year" => first_rank["simulation_end_year"],
            "output_path" => abspath(final_output),
            "partition_count" => count,
            "created_at" => string(now()),
        ))
        # Rank products are diagnostics for a failed merge, not final products.
        # Remove them only once both final products exist.
        isfile(final_output) && isfile(final_manifest) ||
            error("production MPI final products were not written")
        if Bool(get(get(config, "mpi", Dict{String, Any}()), "cleanup_rank_outputs", false))
            rm(joinpath(output_root, ".mpi_ranks", patch); recursive = true, force = true)
        end
    end
    MPI.Barrier(comm)
    return final_manifest
end

if abspath(PROGRAM_FILE) == @__FILE__
    MPI.Init()
    try
        result = run_global_production_mpi()
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 && println(result)
    catch exception
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        println(stderr, "production MPI rank $rank failed: ",
                sprint(showerror, exception, catch_backtrace()))
        MPI.Abort(MPI.COMM_WORLD, 1)
    finally
        MPI.Finalize()
    end
end
