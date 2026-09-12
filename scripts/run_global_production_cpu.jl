# Production-only global entry: run the simulated years from a soil-pool
# allocation that a previous 600-year target-constrained calibration produced.
#
#   julia --project=. scripts/run_global_production_cpu.jl CONFIG.toml CFT_ID rainfed|irrigated
#
# WHY A SEPARATE ENTRY. `run_global_cfts` derives a calibration phase and a
# production phase from one config and always warms up before the production
# years - `agricultural_warmup!` is called unconditionally, and with no warm-up
# keys in the config the year count falls back to ten. A production run that
# continues from a completed calibration has nothing to spin up: `model_inputs`
# applies the calibrated pool allocation to the initial soil state, so those ten
# free years only walk the state away from the calibration it was just handed.
#
# This mirrors the reference workflow, which solves the same problem the same
# way: `work/diff_task/paper1_global_cell/scripts/global_simulation/runtime/
# run_global_historical_final14.jl` reads `[paths] allocation_root` and never
# calls a warm-up. This file is that entry without the reference's own parameter
# machinery, which this project does not use.
#
# It reuses the two-phase runner's helpers rather than copying them - selection,
# management, `model_inputs`, `create_simulation`, `run_output_block!`,
# `write_reconstructed_output`, `partition_cell_selection` - so the only thing
# that differs from a two-phase production run is the absent warm-up.

# The cfts runner brings the wheat runner with it, plus `_file_fingerprint` and
# `_repository_commit`, which the manifest needs for provenance. Guarded so the
# MPI wrapper can include this file after including either of them.
isdefined(@__MODULE__, :_file_fingerprint) ||
    include(joinpath(@__DIR__, "run_global_cfts_cpu.jl"))

"""Resolve the calibration's soil-pool allocation for one patch under `root`.

Copied from the reference entry, including the fallback for MPI products written
before merged patch outputs adopted the `batches/<patch>` layout. Returns the
first direct candidate when nothing exists, so the caller's `isfile` check
reports a path a human can look for.
"""
function allocation_path(root, cft_id, irrigated)
    water = irrigated ? "irrigated" : "rainfed"
    batch = "cft_$(lpad(cft_id, 2, '0'))_$(water)"
    filename = "warmup_soil_pool_allocation.nc"
    direct = [
        joinpath(root, batch, "calibration", filename),
        joinpath(root, "batches", batch, "calibration", filename),
    ]
    for path in direct
        isfile(path) && return path
    end
    nested = String[]
    for batch_root in unique((joinpath(root, batch), joinpath(root, "batches", batch)))
        isdir(batch_root) || continue
        for directory in readdir(batch_root; join = true)
            isdir(directory) && startswith(basename(directory), "mpi_") || continue
            candidate = joinpath(
                directory, "merged", "batches", batch, "calibration", filename,
            )
            isfile(candidate) && push!(nested, candidate)
        end
    end
    length(nested) <= 1 || error("multiple MPI soil allocations found for $batch")
    return isempty(nested) ? first(direct) : only(nested)
end

"""Narrow `selection` to the cells the allocation covers, per the configured
policy.

`strict` refuses any gap: with management fixed at one year the calibration saw
every cell production will, so a gap is a real error. `intersect` drops the
uncovered cells and records how many, which is what transient land use needs -
it can select a cell the calibration never grew a crop in.
"""
function apply_allocation_policy(grid, selection, allocation, policy::Symbol)
    policy in (:strict, :intersect) ||
        error("run.pool_allocation_cell_policy must be strict or intersect")
    covered_ids = Set(allocation.selection.cell_ids)
    covered = map(id -> id in covered_ids, selection.cell_ids)
    excluded = count(!, covered)
    if excluded == 0
        return (selection = selection, excluded = 0)
    end
    policy === :strict && error(
        "pool allocation does not cover every selected cell: $excluded of " *
        "$(length(selection.cell_ids)) are missing. Either the calibration ran on a " *
        "different cell set, or run.pool_allocation_cell_policy should be intersect.",
    )
    narrowed = select_cells(grid, selection.compact_indices[covered])
    isempty(narrowed.cell_ids) &&
        error("pool-allocation intersection removed every selected cell")
    return (selection = narrowed, excluded = excluded)
end

"""One line naming a bad final soil state, with up to eight examples.

Copied from the reference production entry. The examples carry `cell_id`
because a bare index says nothing about which cell to go and look at.
"""
function _soil_state_failure(kind, name, values_host, indices, cell_ids)
    examples = String[]
    for index in Iterators.take(indices, 8)
        coordinates = Tuple(index)
        cell_position = last(coordinates)
        cell_id = if isnothing(cell_ids) || cell_position > length(cell_ids)
            "unknown"
        else
            string(cell_ids[cell_position])
        end
        push!(examples, "index=$coordinates cell_id=$cell_id value=$(values_host[index])")
    end
    return "$kind final soil state in $name; count=$(length(indices)); " *
        "examples=$(join(examples, "; "))"
end

"""Refuse a run whose final soil pools are non-finite or negative.

Copied from the reference production entry. The two-phase runner reaches the
same question through `balance_report` on the diagnostic cells; a production run
has no diagnostic phase, so the check is on the pools themselves and covers
every cell rather than the first ten.
"""
function validate_final_state(simulation; cell_ids = nothing)
    soil = simulation.state.prognostic.soil
    pools = Pair{String, Any}["water.storage" => soil.water.storage]
    append!(pools, ["carbon.$name" => pool for (name, pool) in pairs(soil.carbon)])
    append!(pools, ["nitrogen.$name" => pool for (name, pool) in pairs(soil.nitrogen)])
    for (name, pool) in pools
        values_host = Array(pool)
        if !all(isfinite, values_host)
            indices = findall(value -> !isfinite(value), values_host)
            error(_soil_state_failure("non-finite", name, values_host, indices, cell_ids))
        end
        if !all(>=(0), values_host)
            indices = findall(<(0), values_host)
            error(_soil_state_failure("negative", name, values_host, indices, cell_ids))
        end
    end
end

function run_global_production(
    config_path;
    backend_override = :cpu,
    cft_id::Integer = 1,
    irrigated::Bool = false,
    output_directory_override::Union{Nothing, AbstractString} = nothing,
    partition_rank::Integer = 0,
    partition_count::Integer = 1,
)
    config = TOML.parsefile(config_path)
    run = config["run"]
    paths = config["paths"]
    backend = execution_backend(config; override = backend_override)
    println("execution backend: ", backend.name)

    haskey(paths, "allocation_root") || error(
        "[paths] allocation_root is required: this entry runs production from a " *
        "completed calibration and does not calibrate",
    )
    for key in ("warmup_minimum_years", "warmup_maximum_years", "warmup_years")
        haskey(run, key) && error(
            "[run] $key has no meaning here: this entry never warms up. Remove it, " *
            "or use the two-phase entry.",
        )
    end
    Bool(get(run, "calibration_only", false)) &&
        error("run.calibration_only = true has no meaning for the production entry")

    output_directory = abspath(
        isnothing(output_directory_override) ? paths["output_directory"] :
        output_directory_override,
    )
    mkpath(output_directory)

    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    simulation_years = collect(configured_simulation_years(config))
    start_year = first(simulation_years)
    end_year = last(simulation_years)
    source_years = management_source_years(config, simulation_years)
    nitrogen_deposition_year = nitrogen_deposition_source_year(config)

    # Eligibility, identical to the two-phase entry: land-use fraction, arable
    # soil, and a usable PHU in every active year.
    landuse = read_management(
        catalog, :landuse, grid, cft_id;
        simulation_years = source_years, T = Float32, irrigated,
    )
    crop_mask = build_crop_mask(grid, landuse.values)
    soil_probe = read_soil_data(catalog, grid; selection = crop_mask.selection)
    rock_ice = .!is_arable_crop_soil.(soil_probe.soilcode)
    phu_probe = read_management(
        catalog, :phu, grid, cft_id;
        selection = crop_mask.selection,
        simulation_years = source_years,
        active = falses(size(crop_mask.active)),
        T = Float32, irrigated,
    )
    valid_phu = map(axes(phu_probe.values, 2)) do cell
        active_years = view(crop_mask.active, :, cell)
        values = view(phu_probe.values, :, cell)
        all(!active || (isfinite(value) && value != 0)
            for (active, value) in zip(active_years, values))
    end
    eligible = valid_phu .& .!rock_ice
    any(eligible) || error("no landfrac-selected cells have valid PHU on arable soil")
    selection = select_cells(grid, crop_mask.selection.compact_indices[eligible])

    cell_limit = Int(get(run, "cell_limit", 0))
    cell_limit >= 0 || error("cell_limit must be non-negative")
    if cell_limit > 0
        selection = select_cells(
            grid,
            selection.compact_indices[1:min(cell_limit, length(selection.cell_ids))],
        )
    end
    global_cell_count = length(selection.cell_ids)

    # Rank slicing, through the same helper the two-phase entry uses, so a
    # 32-rank production partitions exactly as a 32-rank calibration did.
    if partition_count > 1
        selection = partition_cell_selection(selection, partition_rank, partition_count)
    end

    allocation_root = abspath(paths["allocation_root"])
    allocation_file = allocation_path(allocation_root, cft_id, irrigated)
    isfile(allocation_file) ||
        error("missing 600-year soil allocation: $allocation_file")
    allocation = read_soil_pool_allocation(allocation_file)
    allocation.cft_id == cft_id || error(
        "allocation CFT metadata mismatch: file says $(allocation.cft_id), run wants $cft_id",
    )
    allocation.irrigated == irrigated || error(
        "allocation irrigation metadata mismatch: file says $(allocation.irrigated)",
    )
    policy = Symbol(lowercase(String(get(run, "pool_allocation_cell_policy", "strict"))))
    narrowed = apply_allocation_policy(grid, selection, allocation, policy)
    selection = narrowed.selection
    excluded_cell_count = narrowed.excluded

    println("production rank $partition_rank/$partition_count: ",
            length(selection.cell_ids), " of $global_cell_count cells, ",
            "allocation=", allocation_file)

    hwsd = load_hwsd_targets(paths, selection)
    selected_landuse = read_management(
        catalog, :landuse, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated,
    )
    management = read_management_schedule(
        catalog, grid, selection, selected_landuse.values .> 0,
        source_years, simulation_years, config, cft_id; irrigated,
    )
    management_blocks = [
        annual_management(management, index) for index in eachindex(simulation_years)
    ]
    initial_data = model_inputs(
        grid, selection, hwsd, catalog, management;
        pool_allocation = allocation,
    )

    expected_days = 365 * length(simulation_years)
    reader = climate_blocks(
        catalog, grid;
        selection, start_year, end_year,
        block_days = 365, T = Float32, nitrogen_deposition_year,
    )
    climate_days(reader) == expected_days ||
        error("expected exactly $expected_days forcing days")
    forcings = climate_forcings(reader)
    length(forcings) == length(simulation_years) ||
        error("365-day blocks must produce one forcing block per simulation year")

    simulation = create_simulation(
        initial_data, selection, config, expected_days, backend.device, cft_id;
        irrigated, diagnostics = false,
    )

    annual_chunks = OutputChunk[]
    compact_writer = NetCDFBlockWriter(joinpath(output_directory, "compact"); prefix = "wheat")
    writer = chunk -> begin
        chunk.frequency === :annual && push!(annual_chunks, chunk)
        compact_writer(chunk)
    end
    stream = OutputStream(
        production_output_variables();
        frequency = :annual, writer, cell_ids = selection.cell_ids,
    )
    memory = estimate_memory(simulation, reader; prefetch = false, output_stream = stream)
    write_report(joinpath(output_directory, "memory_estimate.toml"), Dict(
        string(name) => value for (name, value) in pairs(memory)
    ))

    for index in eachindex(simulation_years)
        run_output_block!(simulation, forcings[index], management_blocks[index], stream)
        simulation.simulated_days == 365 * index || error(
            "$(simulation_years[index]) run did not stop at day $(365 * index)",
        )
    end

    validate_final_state(simulation; cell_ids = selection.cell_ids)

    production_output = joinpath(
        output_directory, "global_wheat_$(start_year)_$(end_year).nc",
    )
    write_reconstructed_output(
        production_output, grid, selection, annual_chunks, simulation_years,
    )

    # The keys the MPI merge needs. `partition_cell_ids` is what lets the merge
    # place each rank's cells and detect an overlap between ranks; without it a
    # merge can only trust the slicing it was told about. `output_path` and the
    # consistency keys are what the reference merge compares across ranks.
    processes = get(config, "processes", Dict{String, Any}())
    management = get(config, "management", Dict{String, Any}())
    manifest = Dict{String, Any}(
        "schema_version" => "1",
        "entry" => "run_global_production_cpu.jl",
        "warmup" => "none: continued from the calibration allocation",
        "repository_commit" => _repository_commit(),
        "config_path" => abspath(config_path),
        "config_fingerprint" => _file_fingerprint(config_path),
        "cft_id" => cft_id,
        "water_system" => irrigated ? "irrigated" : "rainfed",
        "processes_configuration" => String(get(processes, "configuration", "")),
        "rate_scale" => Float64(get(processes, "rate_scale", 1.0)),
        # The RESOLVED value, read back off the CFT the simulation is actually
        # holding, not the config's spec string: "fao56" means a different number
        # for each crop, and a spec cannot show that the value reached the model.
        # `anthesis_heat` once failed exactly here - the keyword was never
        # threaded through the global entry, every unit test passed because they
        # called the kernel directly, and twenty-four jobs ran a configuration
        # nobody had asked for. A manifest that records the number makes that
        # failure visible in the output tree instead of in the yields.
        "depletion_fraction" => Float64(simulation.cft.depletion_fraction),
        "management_mode" => String(get(management, "mode", "")),
        "management_fixed_year" => Int(get(management, "fixed_year", 0)),
        "crop_resp_fix" => Bool(get(run, "crop_resp_fix", false)),
        "nitrogen_limit_vcmax" => Bool(get(run, "nitrogen_limit_vcmax", false)),
        "candidate_cell_count" => global_cell_count,
        "cells" => length(selection.cell_ids),
        "partition_rank" => Int(partition_rank),
        "partition_count" => Int(partition_count),
        "partition_cell_ids" => Int.(selection.cell_ids),
        "simulation_start_year" => start_year,
        "simulation_end_year" => end_year,
        "allocation_root" => allocation_root,
        "allocation_path" => abspath(allocation_file),
        "pool_allocation_cell_policy" => String(policy),
        "allocation_excluded_cell_count" => excluded_cell_count,
        "output_path" => abspath(production_output),
        "backend" => String(Symbol(backend.name)),
    )
    manifest_path = write_report(
        joinpath(output_directory, "run_manifest.toml"), manifest,
    )
    println("output = ", production_output)
    # The manifest is a TOML document, so its keys are strings; the reference
    # returns the two paths rather than the document, and the MPI wrapper
    # re-reads the manifests from disk instead of trusting a return value.
    return (output_path = production_output, manifest_path = manifest_path)
end

"""Positional entry, shaped like the reference's `main` so the MPI wrapper can
hand it a rank, a count and a directory without keyword plumbing."""
function production_main(args = ARGS)
    length(args) in (3, 6) || error(
        "usage: run_global_production_cpu.jl CONFIG_TOML CFT_ID rainfed|irrigated " *
        "[RANK COUNT RANK_DIRECTORY]",
    )
    water = lowercase(args[3])
    water in ("rainfed", "irrigated") ||
        error("water system must be rainfed or irrigated")
    rank = length(args) == 6 ? parse(Int, args[4]) : 0
    count = length(args) == 6 ? parse(Int, args[5]) : 1
    directory = length(args) == 6 ? args[6] : nothing
    return run_global_production(
        abspath(args[1]);
        cft_id = parse(Int, args[2]),
        irrigated = water == "irrigated",
        partition_rank = rank,
        partition_count = count,
        output_directory_override = directory,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(production_main())
end
