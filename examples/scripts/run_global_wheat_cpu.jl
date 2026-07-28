using Agrocosm
using Dates
using NCDatasets
using TOML

include(joinpath(@__DIR__, "..", "..", "lib", "AgrocosmData", "src", "AgrocosmData.jl"))
using .AgrocosmData

function catalog_from_config(config)
    paths = config["paths"]
    input_directory = abspath(paths["input_directory"])
    management_directory = joinpath(input_directory, "management")
    climate_directory = joinpath(input_directory, "climate")
    soil_directory = joinpath(input_directory, "soil")
    registry = PFTRegistry([1], ["temperate cereals"])
    single_pft(filename, variable; units = "") = DatasetSpec(
        joinpath(management_directory, filename), variable; units, pft_ids = [1],
    )
    return DatasetCatalog(
        Dict{Symbol, DatasetSpec}(
            :grid => DatasetSpec(joinpath(soil_directory, "grid.nc"), "cellid"),
            :soilcode => DatasetSpec(joinpath(soil_directory, "soil_30arcmin_13_types.nc"), "soilcode"),
            :soilph => DatasetSpec(joinpath(soil_directory, "soil_pH30arcmin.nc"), "soilph"),
            :landuse => single_pft("landuse_wheat_rainfed.nc", "landfrac"),
            :sowing_date => single_pft("sdate_wheat_rainfed.nc", "sdate"),
            :phu => single_pft("phu_wheat_rainfed.nc", "phusum"),
            :fertilizer => single_pft("fertilizer_wheat_rainfed.nc", "fertilizer"),
            :manure => single_pft("manure_wheat_rainfed.nc", "manure"),
            :residue_fraction => single_pft("residue_wheat_rainfed.nc", "residuefrac"),
            :temp => DatasetSpec(joinpath(climate_directory, "temp_2015_2016.nc"), "temp"),
            :prec => DatasetSpec(joinpath(climate_directory, "prec_2015_2016.nc"), "prec"),
            :lwnet => DatasetSpec(joinpath(climate_directory, "lwnet_2015_2016.nc"), "lwnet"),
            :swdown => DatasetSpec(joinpath(climate_directory, "swdown_2015_2016.nc"), "swdown"),
            :co2 => DatasetSpec(joinpath(climate_directory, "co2_2015_2016.txt"), "co2"),
        ),
        registry,
    )
end

function selected_targets(targets, selection)
    position = Dict(cell_id => index for (index, cell_id) in pairs(targets.selection.cell_ids))
    indices = [get(position, cell_id, 0) for cell_id in selection.cell_ids]
    all(>(0), indices) || error("HWSD targets do not cover every selected wheat cell")
    return SoilCNTargets(
        selection,
        targets.layer_bounds,
        targets.soil_organic_carbon[:, indices],
        targets.total_nitrogen[:, indices],
        targets.coverage[:, indices],
        targets.uncertain[:, indices],
        targets.provenance,
    )
end

function load_hwsd_targets(paths, selection)
    if haskey(paths, "hwsd_targets") && isfile(abspath(paths["hwsd_targets"]))
        return read_soil_cn_targets(abspath(paths["hwsd_targets"]); T = Float32)
    end
    haskey(paths, "hwsd_profile_directory") || error(
        "provide an existing hwsd_targets file or hwsd_profile_directory",
    )
    directory = abspath(paths["hwsd_profile_directory"])
    profiles = [
        read_soil_cn_targets(joinpath(directory, "cell_$(cell_id).nc"); T = Float32)
        for cell_id in selection.cell_ids
    ]
    return SoilCNTargets(
        selection,
        first(profiles).layer_bounds,
        reduce(hcat, (profile.soil_organic_carbon for profile in profiles)),
        reduce(hcat, (profile.total_nitrogen for profile in profiles)),
        reduce(hcat, (profile.coverage for profile in profiles)),
        BitMatrix(reduce(hcat, (profile.uncertain for profile in profiles))),
        (
            schema_version = DATA_SCHEMA_VERSION,
            preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
            source_version = "HWSD v2.01 per-cell profiles",
            deep_rule = :extend_deepest_density,
        ),
    )
end

function model_inputs(catalog, grid, selection, hwsd_targets, config)
    active = trues(1, length(selection.cell_ids))
    sowing_date = read_management(
        catalog, :sowing_date, grid, 1; selection, active, T = Float32,
    )
    phu = read_management(catalog, :phu, grid, 1; selection, active, T = Float32)
    fertilizer = read_management(
        catalog, :fertilizer, grid, 1; selection, T = Float32,
    )
    manure = read_management(catalog, :manure, grid, 1; selection, T = Float32)
    residue = read_management(
        catalog, :residue_fraction, grid, 1; selection, T = Float32,
    )
    management = config["management"]
    crop = crop_inputs(
        ; sowing_date, phu, fertilizer, manure, residue_fraction = residue,
        fertilizer_mode = Symbol(management["fertilizer"]),
        manure_enabled = management["manure"],
    )
    soil = read_soil_data(catalog, grid; selection)
    targets = selected_targets(hwsd_targets, selection)
    initial_state = hwsd_initial_state(targets, soil)
    return model_initial_data(grid, soil, crop, initial_state)
end

function create_simulation(initial_data, selection, config; diagnostics)
    management = config["management"]
    return initialize_simulation(
        cft1, initial_data;
        days = 730,
        indices = collect(1:length(selection.cell_ids)),
        T = Float32,
        diagnostics,
        irrigation = false,
        manure = management["manure"],
        fertilizer = Symbol(management["fertilizer"]),
        with_tillage = management["with_tillage"],
    )
end

function state_equal(left, right)
    left isa AbstractArray && return Array(left) == Array(right)
    typeof(left) == typeof(right) || return false
    fields = fieldnames(typeof(left))
    isempty(fields) && return isequal(left, right)
    return all(name -> state_equal(getfield(left, name), getfield(right, name)), fields)
end

function run_output_block!(simulation, forcing, stream)
    first_day = simulation.simulated_days + 1
    selected = Set((variable.group, variable.field) for variable in stream.variables)
    run_simulation!(
        simulation, forcing;
        spinup = false, reuse_output = true, selected_output = selected,
    )
    rows = simulation.simulated_days - first_day + 1
    consume_output!(stream, simulation.output, first_day; rows)
    clear_output_timeseries!(simulation.output)
    return simulation
end

function write_reconstructed_output(path, grid, selection, chunks)
    years = Int32[2015, 2016]
    length(chunks) == 2 || error("expected two annual output chunks, got $(length(chunks))")
    names = sort!(collect(keys(chunks[1].values)); by = string)
    isfile(path) && rm(path; force = true)
    NCDataset(path, "c") do dataset
        defDim(dataset, "longitude", length(grid.longitude))
        defDim(dataset, "latitude", length(grid.latitude))
        defDim(dataset, "time", 2)
        defVar(dataset, "longitude", Float32, ("longitude",))[:] = grid.longitude
        defVar(dataset, "latitude", Float32, ("latitude",))[:] = grid.latitude
        defVar(dataset, "time", Int32, ("time",))[:] = years
        defVar(dataset, "cellid", Int32, ("longitude", "latitude"))[:, :] = grid.cellid
        for name in names
            values = Array{Float32}(undef, length(grid.longitude), length(grid.latitude), 2)
            for year_index in 1:2
                values[:, :, year_index] = expand_to_grid(
                    vec(chunks[year_index].values[name]), grid; selection,
                )
            end
            variable = defVar(dataset, String(name), Float32, ("longitude", "latitude", "time"))
            variable[:, :, :] = values
        end
    end
    return path
end

toml_value(value::NamedTuple) = Dict(string(name) => toml_value(item) for (name, item) in pairs(value))
toml_value(value::AbstractDict) = Dict(string(name) => toml_value(item) for (name, item) in pairs(value))
toml_value(value::AbstractArray) = toml_value.(collect(value))
toml_value(value::DataType) = string(value)
toml_value(value::Symbol) = string(value)
toml_value(value) = value

function write_report(path, report)
    open(path, "w") do io
        TOML.print(io, toml_value(report); sorted = true)
    end
    return path
end

function balance_report(simulation, validation)
    water = maximum(abs, Array(simulation.water_balance.residual))
    carbon = maximum(abs, Array(simulation.carbon_balance.relative_residual))
    nitrogen = maximum(abs, Array(simulation.nitrogen_balance.relative_residual))
    energy = Array(simulation.thermal_balance.energy_residual)
    column_energy = Array(simulation.thermal_balance.column_energy)
    thermal = maximum(abs.(energy) ./ max.(abs.(column_energy), 1.0f0))
    percolation = maximum(abs, Array(simulation.thermal_balance.percolation_energy_residual))
    water <= validation["maximum_water_residual"] || error("sampled water closure failed")
    carbon <= validation["maximum_carbon_relative_residual"] || error("sampled carbon closure failed")
    nitrogen <= validation["maximum_nitrogen_relative_residual"] || error("sampled nitrogen closure failed")
    thermal <= validation["maximum_thermal_relative_residual"] || error("sampled thermal closure failed")
    percolation <= validation["maximum_percolation_energy_residual"] ||
        error("sampled percolation-energy closure failed")
    return Dict{String, Any}(
        "maximum_water_residual" => water,
        "maximum_carbon_relative_residual" => carbon,
        "maximum_nitrogen_relative_residual" => nitrogen,
        "maximum_thermal_relative_residual" => thermal,
        "maximum_percolation_energy_residual" => percolation,
    )
end

function run_global_wheat(config_path)
    config = TOML.parsefile(config_path)
    paths = config["paths"]
    run = config["run"]
    output_directory = abspath(paths["output_directory"])
    mkpath(output_directory)
    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    landuse = read_management(catalog, :landuse, grid, 1; T = Float32)
    size(landuse.values, 1) == 1 || error("landuse must contain only fixed 2015 management")
    selection = build_crop_mask(grid, landuse.values).selection
    cell_limit = Int(get(run, "cell_limit", 0))
    cell_limit >= 0 || error("cell_limit must be non-negative")
    if cell_limit > 0
        selection = select_cells(
            grid, selection.compact_indices[1:min(cell_limit, length(selection.cell_ids))],
        )
    end
    hwsd = load_hwsd_targets(paths, selection)
    initial_data = model_inputs(catalog, grid, selection, hwsd, config)

    reader = climate_blocks(catalog, grid; selection, block_days = 365, T = Float32)
    climate_days(reader) == 730 || error("expected exactly 730 forcing days")
    first_block = read_climate_block(reader, 1)
    last_block = read_climate_block(reader, length(reader))
    start_date = first(first_block.time)
    end_date = last(last_block.time)
    (year(start_date), month(start_date), day(start_date)) == (2015, 1, 1) ||
        error("forcing must start on 2015-01-01, got $start_date")
    (year(end_date), month(end_date), day(end_date)) == (2016, 12, 31) ||
        error("forcing must end on 2016-12-31, got $end_date")
    forcings = climate_forcings(reader)
    length(forcings) == 2 || error("365-day blocks must produce exactly two forcing blocks")

    warmup_simulation = create_simulation(initial_data, selection, config; diagnostics = false)
    warmup = agricultural_warmup!(
        warmup_simulation, forcings; years = Int(run["warmup_years"]),
    )
    warmup_checkpoint = joinpath(output_directory, "warmup_checkpoint.jld2")
    save_checkpoint(warmup_checkpoint, warmup_simulation)
    write_report(
        joinpath(output_directory, "warmup_cn_drift.toml"),
        agricultural_warmup_drift(warmup),
    )

    production = create_simulation(initial_data, selection, config; diagnostics = false)
    restore_checkpoint!(production, warmup_checkpoint)
    state_equal(production.state.prognostic, warmup_simulation.state.prognostic) ||
        error("warm-up checkpoint did not restore the prognostic state exactly")

    annual_chunks = OutputChunk[]
    compact_writer = NetCDFBlockWriter(joinpath(output_directory, "compact"); prefix = "wheat")
    writer = chunk -> begin
        chunk.frequency === :annual && push!(annual_chunks, chunk)
        compact_writer(chunk)
    end
    stream = OutputStream(
        [
            OutputVariable(:crop, :gpp; reduction = :sum),
            OutputVariable(:crop, :npp; reduction = :sum),
            OutputVariable(:crop, :lai; reduction = :mean),
            OutputVariable(:crop, :yield),
        ];
        frequency = :annual, writer, cell_ids = selection.cell_ids,
    )
    memory = estimate_memory(
        production, reader; prefetch = false, output_stream = stream,
    )
    write_report(joinpath(output_directory, "memory_estimate.toml"), Dict(
        string(name) => value for (name, value) in pairs(memory)
    ))

    run_output_block!(production, forcings[1], stream)
    production.simulated_days == 365 || error("2015 run did not stop at day 365")
    year2015_checkpoint = joinpath(output_directory, "end_2015_checkpoint.jld2")
    save_checkpoint(year2015_checkpoint, production)

    continued = create_simulation(initial_data, selection, config; diagnostics = false)
    restore_checkpoint!(continued, year2015_checkpoint)
    state_equal(continued.state.prognostic, production.state.prognostic) ||
        error("2015 checkpoint did not restore the prognostic state exactly")
    run_output_block!(continued, forcings[2], stream)
    finish_output_stream!(stream, continued.simulated_days)
    continued.simulated_days == 730 || error("production run did not reach the end of 2016")
    all(isfinite, continued.state.prognostic.soil.water.storage) || error("non-finite soil water")
    all(>=(0), continued.state.prognostic.soil.water.storage) || error("negative soil water")
    for pool in values(continued.state.prognostic.soil.carbon)
        all(isfinite, pool) && all(>=(0), pool) || error("invalid soil carbon pool")
    end
    for pool in values(continued.state.prognostic.soil.nitrogen)
        all(isfinite, pool) && all(>=(0), pool) || error("invalid soil nitrogen pool")
    end
    write_reconstructed_output(
        joinpath(output_directory, "global_wheat_2015_2016.nc"),
        grid, selection, annual_chunks,
    )

    diagnostic_count = min(Int(run["diagnostic_cells"]), length(selection.cell_ids))
    diagnostic_selection = select_cells(grid, selection.compact_indices[1:diagnostic_count])
    diagnostic_initial = model_inputs(catalog, grid, diagnostic_selection, hwsd, config)
    diagnostic_reader = climate_blocks(
        catalog, grid; selection = diagnostic_selection, block_days = 365, T = Float32,
    )
    diagnostic = create_simulation(
        diagnostic_initial, diagnostic_selection, config; diagnostics = true,
    )
    agricultural_warmup!(
        diagnostic, climate_forcings(diagnostic_reader); years = Int(run["warmup_years"]),
    )
    run_simulation!(diagnostic, climate_forcings(diagnostic_reader); spinup = false)
    balance = balance_report(diagnostic, config["validation"])
    balance["simulation_summary"] = simulation_summary(diagnostic)
    write_report(joinpath(output_directory, "sampled_balance_summary.toml"), balance)
    return (; cells = length(selection.cell_ids), memory, output_directory)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_wheat_cpu.jl CONFIG_TOML")
    println(run_global_wheat(abspath(ARGS[1])))
end
