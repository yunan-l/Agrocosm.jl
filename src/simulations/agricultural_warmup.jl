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

    snapshots = Vector{typeof(_warmup_soil_snapshot(simulation.state))}(undef, years)
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
        soil = _warmup_history(snapshots),
    )
end
