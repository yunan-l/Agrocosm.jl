# Reuse the model's input preparation and output contracts, not a trained
# parameter-selection workflow. This file never launches a simulation itself.
include(joinpath(@__DIR__, "..", "run_global_wheat_cpu.jl"))

function write_compact_history(path, grid, selection, chunks, years)
    ispath(path) && error("refusing to overwrite $path")
    length(chunks) == length(years) || error("annual output length mismatch")
    longitude_indices = grid.longitude_indices[selection.compact_indices]
    latitude_indices = grid.latitude_indices[selection.compact_indices]
    NCDataset(path, "c") do dataset
        dataset.attrib["calendar"] = "365_day"
        dataset.attrib["parameter_source"] = "version_defaults"
        defDim(dataset, "cell", length(selection.cell_ids))
        defDim(dataset, "year", length(years))
        defVar(dataset, "cell_id", Int32, ("cell",))[:] = selection.cell_ids
        defVar(dataset, "year", Int32, ("year",))[:] = years
        defVar(dataset, "longitude", Float32, ("cell",))[:] =
            grid.longitude[longitude_indices]
        defVar(dataset, "latitude", Float32, ("cell",))[:] =
            grid.latitude[latitude_indices]
        for variable in production_output_variables()
            name = Symbol(variable.group, :_, variable.field)
            spec, _ = Agrocosm.output_variable_spec(variable.group, variable.field)
            output = defVar(dataset, String(name), Float32, ("cell", "year");
                attrib = Dict("units" => spec.units))
            for (index, chunk) in enumerate(chunks)
                values = vec(chunk.values[name])
                all(isfinite, values) || error("non-finite annual $name in year $(years[index])")
                name == :crop_yield && any(<(0), values) && error("negative harvest yield")
                output[:, index] = values
            end
        end
    end
    return path
end

function default_historical_context(config, config_path, cft_id, water;
    rank = 0, ranks = 1, cell_ids = nothing, end_year = 2019,
)
    record = read_stage_completion(config, config_path, "warmup_600y", cft_id, water)
    allocation_file = require_inside(record["allocation_path"], config["study"]["results_root"])
    file_sha(allocation_file) == record["allocation_sha256"] || error("warm-up allocation has changed")
    allocation = read_soil_pool_allocation(allocation_file)
    irrigated = water == "irrigated"
    allocation.cft_id == cft_id && allocation.irrigated == irrigated || error("warm-up patch mismatch")
    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    covered = Set(allocation.selection.cell_ids)
    wanted = cell_ids === nothing ? covered : Set(Int32.(cell_ids))
    issubset(wanted, covered) || error("requested cells are not covered by the new warm-up")
    selection = select_cells(grid, map(id -> id in wanted, grid.cell_ids))
    length(selection.cell_ids) == length(wanted) || error("canonical grid does not cover requested cells")
    ranks == 1 || (selection = partition_cell_selection(selection, rank, ranks))
    years = collect(1901:end_year)
    source_years = management_source_years(config, years)
    landuse = read_management(catalog, :landuse, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated)
    management = read_management_schedule(catalog, grid, selection, landuse.values .> 0,
        source_years, years, config, cft_id; irrigated)
    management_blocks = [annual_management(management, index) for index in eachindex(years)]
    initial = model_inputs(grid, selection, load_hwsd_targets(config["paths"], selection),
        catalog, management; pool_allocation = allocation)
    reader = climate_blocks(catalog, grid; selection, start_year = 1901, end_year,
        block_days = 365, T = Float32,
        nitrogen_deposition_year = nitrogen_deposition_source_year(config))
    climate_days(reader) == 365 * length(years) || error("historical forcing is incomplete")
    simulation = create_simulation(initial, selection, config, 365 * length(years),
        identity, cft_id; irrigated, diagnostics = false)
    return (; simulation, grid, selection, years, management_blocks,
        forcings = climate_forcings(reader), allocation_file)
end
