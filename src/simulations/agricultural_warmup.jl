function _warmup_soil_snapshot(state::ModelState)
    carbon = soil_carbon_prognostic(state)
    nitrogen = soil_nitrogen_prognostic(state)
    water = soil_water_prognostic(state)
    host(values) = Array(values)

    litter_carbon = vec(sum(host(carbon.litter); dims = 1))
    litter_nitrogen = vec(sum(host(nitrogen.litter); dims = 1))
    fast_carbon = vec(sum(host(carbon.fast); dims = 1))
    slow_carbon = vec(sum(host(carbon.slow); dims = 1))
    fast_nitrogen = vec(sum(host(nitrogen.fast); dims = 1))
    slow_nitrogen = vec(sum(host(nitrogen.slow); dims = 1))
    mineral_nitrogen = vec(sum(
        host(nitrogen.nitrate) .+ host(nitrogen.ammonium); dims = 1,
    ))
    return (
        total_carbon = litter_carbon .+ fast_carbon .+ slow_carbon,
        total_nitrogen = litter_nitrogen .+ fast_nitrogen .+
            slow_nitrogen .+ mineral_nitrogen,
        litter_carbon = litter_carbon,
        litter_nitrogen = litter_nitrogen,
        fast_carbon = fast_carbon,
        fast_nitrogen = fast_nitrogen,
        slow_carbon = slow_carbon,
        slow_nitrogen = slow_nitrogen,
        mineral_nitrogen = mineral_nitrogen,
        soil_water = vec(sum(host(water.storage); dims = 1)),
    )
end

function _warmup_history(snapshots)
    names = keys(first(snapshots))
    return NamedTuple{names}(map(
        name -> reduce(vcat, (permutedims(getproperty(snapshot, name)) for snapshot in snapshots)),
        names,
    ))
end

"""
    agricultural_warmup!(simulation, climate; years=10)

Cycle complete agricultural forcing years before the production run. The
forcing must contain one or more whole 365-day years; multi-year blocks cycle
in order when `years` exceeds the available history. Crop, soil, water,
thermal, carbon, and nitrogen processes are active, but
production outputs, balance ledgers, and `simulation.simulated_days` are left
untouched. The warmed prognostic state is retained.

The returned report contains one host-side row per warm-up year and one column
per active cell. This is a finite transient warm-up, not an equilibrium soil
organic-carbon spin-up.
"""
function agricultural_warmup!(
    simulation::CropSimulation,
    climate::NamedTuple;
    years::Integer = 10,
)
    years > 0 || throw(ArgumentError("warm-up years must be positive"))
    simulation.simulated_days == 0 || throw(ArgumentError(
        "agricultural warm-up must run before the production simulation",
    ))
    _output_timeseries_empty(simulation.output) || throw(ArgumentError(
        "agricultural warm-up requires empty production output",
    ))

    prepared_climate = _prepare_climate(simulation, climate)
    climate_days = size(prepared_climate.temp, 1)
    climate_days > 0 && climate_days % 365 == 0 || throw(ArgumentError(
        "agricultural warm-up requires complete 365-day climate years, got $climate_days rows",
    ))
    cells = length(simulation.config.execution.domain.indices)
    size(prepared_climate.temp, 2) == cells || throw(DimensionMismatch(
        "warm-up climate has $(size(prepared_climate.temp, 2)) cells; simulation has $cells",
    ))

    initial_soil = _warmup_soil_snapshot(simulation.state)
    snapshots = Vector{typeof(initial_soil)}(undef, years)
    no_output = Set{Tuple{Symbol, Symbol}}()
    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        fertilizer = simulation.config.fertilizer,
        with_tillage = simulation.config.with_tillage,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        reuse_output = true,
        selected_output = no_output,
    )

    forcing_years = div(climate_days, 365)
    for year in 1:years
        forcing_year = mod(year - 1, forcing_years)
        start_day = 365 * forcing_year + 1
        end_day = start_day + 364
        _daily_crop!(
            pathway_value(simulation.processes.pathway),
            start_day, end_day, simulation.processes, prepared_climate, simulation.state;
            simulation_day_offset = 365 * (year - 1) + 1 - start_day,
            common...,
        )
        snapshots[year] = _warmup_soil_snapshot(simulation.state)
    end

    # `update_climbuf!` normally closes a year on the following day 1. The
    # production clock intentionally restarts at day 1, so close the final
    # warm-up year explicitly before returning.
    annual_climbuf!(simulation.climbuf.atemp, simulation.climbuf, simulation.pft)
    clear_output_timeseries!(simulation.output)
    return (
        years = Int(years),
        days = 365 * Int(years),
        forcing_years = forcing_years,
        initial_soil = initial_soil,
        soil = _warmup_history(snapshots),
    )
end

"""
    agricultural_warmup!(simulation, climate_blocks; years=10)

Run the agricultural warm-up from restartable climate blocks without joining
the forcing into one in-memory array. Block boundaries may occur anywhere
within a 365-day forcing year.
"""
function agricultural_warmup!(
    simulation::CropSimulation,
    climate_blocks::AbstractVector;
    years::Integer = 10,
)
    years > 0 || throw(ArgumentError("warm-up years must be positive"))
    simulation.simulated_days == 0 || throw(ArgumentError(
        "agricultural warm-up must run before the production simulation",
    ))
    _output_timeseries_empty(simulation.output) || throw(ArgumentError(
        "agricultural warm-up requires empty production output",
    ))
    isempty(climate_blocks) && throw(ArgumentError(
        "agricultural warm-up requires at least one climate block",
    ))

    block_days = map(eachindex(climate_blocks)) do index
        block = climate_blocks[index]
        hasproperty(block, :temp) || throw(ArgumentError(
            "warm-up climate block $index has no temp field",
        ))
        size(block.temp, 1)
    end
    climate_days = sum(block_days)
    climate_days > 0 && climate_days % 365 == 0 || throw(ArgumentError(
        "agricultural warm-up requires complete 365-day climate years, got $climate_days rows",
    ))

    initial_soil = _warmup_soil_snapshot(simulation.state)
    snapshots = Vector{typeof(initial_soil)}(undef, years)
    no_output = Set{Tuple{Symbol, Symbol}}()
    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        fertilizer = simulation.config.fertilizer,
        with_tillage = simulation.config.with_tillage,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        reuse_output = true,
        selected_output = no_output,
    )

    cells = length(simulation.config.execution.domain.indices)
    forcing_years = div(climate_days, 365)
    block_ends = cumsum(block_days)
    warmup_day = 0
    for year in 1:years
        forcing_year = mod(year - 1, forcing_years)
        forcing_start = 365 * forcing_year + 1
        forcing_end = forcing_start + 364

        for index in eachindex(climate_blocks)
            block_start = block_ends[index] - block_days[index] + 1
            block_end = block_ends[index]
            segment_start = max(forcing_start, block_start)
            segment_end = min(forcing_end, block_end)
            segment_start <= segment_end || continue

            prepared_climate = _prepare_climate(simulation, climate_blocks[index])
            size(prepared_climate.temp, 2) == cells || throw(DimensionMismatch(
                "warm-up climate block $index has $(size(prepared_climate.temp, 2)) cells; simulation has $cells",
            ))
            local_start = segment_start - block_start + 1
            local_end = segment_end - block_start + 1
            _daily_crop!(
                pathway_value(simulation.processes.pathway),
                local_start, local_end, simulation.processes, prepared_climate,
                simulation.state;
                simulation_day_offset = warmup_day + 1 - local_start,
                common...,
            )
            warmup_day += local_end - local_start + 1
        end
        snapshots[year] = _warmup_soil_snapshot(simulation.state)
    end

    annual_climbuf!(simulation.climbuf.atemp, simulation.climbuf, simulation.pft)
    clear_output_timeseries!(simulation.output)
    return (
        years = Int(years),
        days = 365 * Int(years),
        forcing_years = forcing_years,
        initial_soil = initial_soil,
        soil = _warmup_history(snapshots),
    )
end
