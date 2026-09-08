# Host-side, selected-cell preparation. None of this runs inside AD or a GPU
# process kernel; the numerical trajectory uses the ordinary model API.
function prepare_historical_case(config, config_path, specification)
    cft_id, water = specification["cft_id"], specification["water_system"]
    year, cell_id = specification["event_year"], specification["cell_id"]
    1902 <= year <= 2019 || error("event year must allow a preceding climate year")
    historical = read_stage_completion(config, config_path, "historical_1901_2019", cft_id, water)
    root = config["study"]["results_root"]
    records = [TOML.parsefile(require_inside(path, root)) for path in historical["rank_manifests"]]
    rank_record = only(filter(record -> cell_id in record["cell_ids"], records))
    rank_record["simulation_config_sha256"] == simulation_config_sha(config) || error("historical rank config mismatch")
    output_path = require_inside(rank_record["output_path"], root)
    file_sha(output_path) == rank_record["output_sha256"] || error("historical output has changed")
    historical_yield = NCDataset(output_path, "r") do dataset
        cell_index = only(findall(==(cell_id), dataset["cell_id"][:]))
        year_index = only(findall(==(year), dataset["year"][:]))
        Float64(dataset["crop_yield"][cell_index, year_index]) / 0.45 * 0.01
    end
    context = default_historical_context(config, config_path, cft_id, water;
        cell_ids = [cell_id], end_year = year)
    simulation = context.simulation
    for index in 1:(length(context.years) - 2)
        run_simulation!(simulation, context.forcings[index]; spinup = false,
            reuse_output = true, management = context.management_blocks[index])
        clear_output_timeseries!(simulation.output)
    end
    previous, current = context.forcings[end - 1], context.forcings[end]
    climate = (
        temp = vcat(previous.temp, current.temp), prec = vcat(previous.prec, current.prec),
        sw = vcat(previous.sw, current.sw), lw = vcat(previous.lw, current.lw),
        wind = vcat(previous.wind, current.wind),
        no3_deposition = vcat(previous.no3_deposition, current.no3_deposition),
        nh4_deposition = vcat(previous.nh4_deposition, current.nh4_deposition),
        co2 = vcat(previous.co2, current.co2), co2_daily = true,
    )
    initial = nothing
    first_day = 0
    season = nothing
    for (offset, forcing) in enumerate((previous, current))
        management = context.management_blocks[end - 2 + offset]
        for day in 1:365
            clear_output_timeseries!(simulation.output)
            transition_day!(simulation, forcing; climate_day = day, management)
            joined_day = 365 * (offset - 1) + day
            if simulation.output.calendar.sowing_event[end, 1] != 0
                initial = deepcopy(simulation.state)
                first_day = joined_day + 1 # The discrete establishment day is held fixed.
            end
            if simulation.output.calendar.harvest_event[end, 1] != 0
                if offset == 2
                    season === nothing || error("multiple target-year harvests require explicit season selection")
                    initial === nothing && error("harvest has no sowing in the two-year replay window")
                    Agrocosm.crop_events(simulation.state).harvest[1] == 0 ||
                        error("target crop failed: fixed-calendar AD cannot describe this event")
                    first_day < joined_day || error("empty established-crop window")
                    season = (; state = initial, days = first_day:(joined_day - 1), harvest_day = joined_day,
                        production_yield = Agrocosm.crop_fluxes(simulation.state).carbon.yield[1] / 0.45f0 * 0.01f0)
                end
                initial = nothing
            end
        end
    end
    season === nothing && error("no target-year harvest; screen another case or analyse the failure separately")
    season.production_yield > 0 || error("zero-yield seasons require a separate crop-failure analysis")
    isapprox(season.production_yield, historical_yield; rtol = 1e-5, atol = 1e-6) ||
        error("selected-cell replay does not reproduce the new historical yield: $(season.production_yield) vs $historical_yield")
    Agrocosm.enzyme_prepare_daily_state!(season.state)
    return merge(season, (; climate, forcing = Agrocosm.weather_forcing(climate),
        cft = simulation.cft, parameters = simulation.model_parameters,
        anchor_year = year - 1, historical_yield, allocation_file = context.allocation_file,
        selection = context.selection, grid = context.grid))
end

function reference_weather(config, case, reference_year)
    1902 <= reference_year <= 2019 || error("reference year must allow a preceding climate year")
    catalog = catalog_from_config(config)
    reader = climate_blocks(catalog, case.grid; selection = case.selection,
        start_year = reference_year - 1, end_year = reference_year, block_days = 365,
        T = eltype(case.forcing), nitrogen_deposition_year = nitrogen_deposition_source_year(config))
    blocks = climate_forcings(reader)
    first_block, second_block = blocks[1], blocks[2]
    controls = vcat(Agrocosm.weather_forcing(first_block), Agrocosm.weather_forcing(second_block))
    result = copy(case.forcing)
    # Replace coherent, calendar-matched weather from post-sowing onward.
    # CO₂/deposition, management and the pre-event state remain factual.
    result[first(case.days):end, :, :] .= controls[first(case.days):end, :, :]
    return result
end

"""Restore coherent reference weather only inside an explicit factual window."""
function weather_window_controls(forcing, reference, days::UnitRange{Int})
    size(forcing) == size(reference) || throw(DimensionMismatch("reference weather shape mismatch"))
    !isempty(days) && first(days) >= 1 && last(days) <= size(forcing, 1) ||
        throw(ArgumentError("invalid weather replacement window"))
    controls = copy(forcing)
    controls[days, :, :] .= reference[days, :, :]
    return controls
end

function write_weather_case_output(path, case, result, factual, controls, counterfactual;
    replacement_days = first(case.days):size(case.forcing, 1),
    counterfactual_kind = "full_reference", reference_year = 0,
)
    ispath(path) && error("refusing to overwrite $path")
    NCDataset(path, "c") do dataset
        dataset.attrib["calendar"] = "365_day"
        dataset.attrib["anchor_year"] = case.anchor_year
        dataset.attrib["cell_id"] = Int(only(case.selection.cell_ids))
        dataset.attrib["yield_units"] = "t dry matter ha-1"
        dataset.attrib["gradient_scope"] = "post-sowing weather; fixed harvest; fixed pre-event state"
        dataset.attrib["process_diagnostics"] = "mediator diagnostics, not additive causal process contributions"
        dataset.attrib["diagnostic_sampling"] = "end-of-day stocks; daily fluxes and last process-stage auxiliaries; exclude harvest reset days"
        dataset.attrib["schema_version"] = 2
        dataset.attrib["counterfactual_kind"] = counterfactual_kind
        dataset.attrib["reference_year"] = reference_year
        dataset.attrib["factual_yield"] = Float64(factual.yield)
        dataset.attrib["counterfactual_yield"] = Float64(counterfactual.yield)
        dataset.attrib["factual_harvest_day"] = case.harvest_day
        dataset.attrib["counterfactual_harvest_day"] = isempty(counterfactual.harvest_days) ? -1 : only(counterfactual.harvest_days)
        dataset.attrib["counterfactual_failed"] = Int(counterfactual.failed)
        dataset.attrib["schedule_matches"] = Int(counterfactual.schedule_matches)
        dataset.attrib["replacement_first_day"] = first(replacement_days)
        dataset.attrib["replacement_last_day"] = last(replacement_days)
        defDim(dataset, "day", size(case.forcing, 1))
        defVar(dataset, "day", Int32, ("day",))[:] = 1:size(case.forcing, 1)
        defVar(dataset, "calendar_year", Int32, ("day",))[:] =
            [case.anchor_year + fld(day - 1, 365) for day in 1:size(case.forcing, 1)]
        defVar(dataset, "day_of_year", Int32, ("day",))[:] = mod1.(1:size(case.forcing, 1), 365)
        active = zeros(Int8, size(case.forcing, 1))
        active[case.days] .= 1
        defVar(dataset, "ad_active_window", Int8, ("day",);
            attrib = Dict("description" => "zero gradients outside this mask are fixed by design"))[:] = active
        replacement = zeros(Int8, size(case.forcing, 1))
        replacement[replacement_days] .= 1
        defVar(dataset, "weather_replacement_window", Int8, ("day",);
            attrib = Dict("description" => "all five weather channels are restored inside this mask only"))[:] = replacement
        for (index, (name, unit)) in enumerate(zip(Agrocosm.WEATHER_VARIABLES,
            ("degree_C", "mm day-1", "W m-2", "W m-2", "m s-1")))
            defVar(dataset, "factual_$(name)", Float32, ("day",); attrib = Dict("units" => unit))[:] = case.forcing[:, 1, index]
            defVar(dataset, "reference_$(name)", Float32, ("day",); attrib = Dict("units" => unit))[:] = controls[:, 1, index]
            defVar(dataset, "dyield_d$(name)", Float32, ("day",);
                attrib = Dict("units" => "(t dry matter ha-1) / ($unit)"))[:] = result.gradient[:, 1, index]
        end
        units = (harvest_event = "1", gpp = "gC m-2 day-1", npp = "gC m-2 day-1",
            respiration = "gC m-2 day-1", biological_fixation_cost = "gC m-2 day-1",
            et = "mm day-1", transpiration = "mm day-1", soil_evaporation = "mm day-1",
            lai = "m2 m-2", apar = "J m-2 day-1", fphu = "1", biomass_carbon = "gC m-2",
            leaf_carbon = "gC m-2", root_carbon = "gC m-2", mobile_carbon = "gC m-2",
            storage_carbon = "gC m-2", soil_water = "mm", rootzone_available_water = "mm",
            water_sufficiency = "1", canopy_conductance = "mm s-1", temperature_stress = "1",
            topsoil_temperature = "degree_C", soil_nitrate = "gN m-2", soil_ammonium = "gN m-2",
            nitrogen_uptake = "gN m-2 day-1", fertilizer_input = "gN m-2 day-1",
            manure_input = "gN m-2 day-1", mineralization = "gN m-2 day-1",
            leaching = "gN m-2 day-1", volatilization = "gN m-2 day-1", leaf_nitrogen = "gN m-2",
            potential_vcmax = "gC m-2 day-1", vcmax = "gC m-2 day-1", nitrogen_limitation = "1")
        for (prefix, replay) in (("factual", factual), ("reference", counterfactual))
            for name in keys(replay.daily)
                name == :day && continue
                values = fill(Float32(NaN), size(case.forcing, 1))
                values[replay.daily.day] .= getproperty(replay.daily, name)
                kind = name in (:gpp, :npp, :respiration, :biological_fixation_cost,
                    :et, :transpiration, :soil_evaporation, :apar, :nitrogen_uptake,
                    :fertilizer_input, :manure_input, :mineralization, :leaching, :volatilization) ? "daily_flux" : "state_or_auxiliary"
                defVar(dataset, "$(prefix)_$(name)", Float32, ("day",);
                    attrib = Dict("units" => getproperty(units, name), "diagnostic_kind" => kind))[:] = values
            end
        end
    end
    return path
end

function write_weather_counterfactual(path, case, result, factual, controls, counterfactual;
    replacement_days, counterfactual_kind, reference_year,
)
    write_weather_case_output(path, case, result, factual, controls, counterfactual;
        replacement_days, counterfactual_kind, reference_year)
    projected_change = sum(Float64.(result.gradient) .* Float64.(controls .- case.forcing))
    return Dict{String, Any}(
        "reference_year" => reference_year, "counterfactual_kind" => counterfactual_kind,
        "replacement_first_day" => first(replacement_days), "replacement_last_day" => last(replacement_days),
        "harvest_observed" => !isempty(counterfactual.harvest_days),
        "yield" => counterfactual.yield, "yield_change" => counterfactual.yield - factual.yield,
        "harvest_days" => counterfactual.harvest_days, "crop_failed" => counterfactual.failed,
        "schedule_matches" => counterfactual.schedule_matches,
        "linearized_fixed_event_yield_change" => projected_change,
        "linearization_residual" => counterfactual.yield - factual.yield - projected_change,
        "output_path" => path, "output_sha256" => file_sha(path),
    )
end
