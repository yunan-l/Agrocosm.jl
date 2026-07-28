using Agrocosm
using JLD2
using Test

include(joinpath(@__DIR__, "..", "..", "lib", "AgrocosmData", "src", "AgrocosmData.jl"))
import .AgrocosmData

function fixture_array_snapshot(value, path = "")
    arrays = Dict{String, Any}()
    if value isa AbstractArray
        arrays[path] = Array(value)
        return arrays
    end
    for name in fieldnames(typeof(value))
        child = getfield(value, name)
        child_path = isempty(path) ? string(name) : string(path, ".", name)
        merge!(arrays, fixture_array_snapshot(child, child_path))
    end
    return arrays
end

function test_fixture_arrays_equal(new_value, old_value)
    new_arrays = fixture_array_snapshot(new_value)
    old_arrays = fixture_array_snapshot(old_value)
    @test keys(new_arrays) == keys(old_arrays)
    for path in keys(old_arrays)
        @test new_arrays[path] == old_arrays[path]
    end
end

@testset "AgrocosmData ten-cell fixture matches legacy loaders" begin
    initial = load(joinpath(@__DIR__, "..", "..", "examples", "initial_wheat.jld2"), "initial_data")
    climate = load(joinpath(@__DIR__, "..", "..", "examples", "climate_2000_2009.jld2"), "climate")
    cells = length(initial.latitude)
    days = 365
    indices = collect(1:cells)

    legacy = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, fertilizer = :yes,
    )
    run_simulation!(legacy, climate; end_day = days, spinup = false)

    selection = AgrocosmData.CellSelection(indices, Int32.(0:(cells - 1)))
    grid = AgrocosmData.GridIndex(
        Float32[0],
        copy(initial.latitude),
        reshape(copy(selection.cell_ids), 1, :),
        copy(selection.cell_ids),
        ones(Int32, cells),
        Int32.(indices),
    )
    soil = AgrocosmData.SoilData(
        selection,
        Int32.(initial.soilparam.soilcode),
        copy(initial.soilparam.soilph),
        copy(initial.soilparam.w_sat),
        copy(initial.soilparam.sand),
        copy(initial.soilparam.silt),
        copy(initial.soilparam.clay),
        copy(initial.soilparam.tdiff_0),
        copy(initial.soilparam.tdiff_15),
        copy(initial.soilparam.soildepth),
        (fixture = "examples/initial_wheat.jld2",),
    )
    state_fields = (:swc, :litc, :fastc, :slowc, :litn, :fastn, :slown)
    initial_state = NamedTuple{state_fields}(map(
        name -> copy(getproperty(initial.initialLPJmL.u0, name)), state_fields,
    ))
    prepared = AgrocosmData.model_initial_data(grid, soil, initial.crop, initial_state)
    adapted = initialize_simulation(
        cft1, prepared;
        T = Float32, days, fertilizer = :yes,
    )

    block_size = 73
    blocks = AgrocosmData.ClimateBlock[]
    for first_day in 1:block_size:days
        last_day = min(days, first_day + block_size - 1)
        rows = first_day:last_day
        push!(blocks, AgrocosmData.ClimateBlock(
            collect(rows),
            copy(climate.temp[rows, :]),
            copy(climate.prec[rows, :]),
            copy(climate.swdown[rows, :]),
            copy(climate.lwnet[rows, :]),
            fill(Float32(climate.co2[1]), length(rows)),
            selection,
            (fixture = "examples/climate_2000_2009.jld2",),
        ))
    end
    run_simulation!(adapted, AgrocosmData.climate_forcings(blocks); spinup = false)

    @test adapted.simulated_days == legacy.simulated_days == days
    test_fixture_arrays_equal(adapted.state.prognostic, legacy.state.prognostic)
    test_fixture_arrays_equal(adapted.state.fluxes, legacy.state.fluxes)
    test_fixture_arrays_equal(adapted.state.auxiliary, legacy.state.auxiliary)
    test_fixture_arrays_equal(adapted.state.events, legacy.state.events)
    test_fixture_arrays_equal(adapted.state.inputs.crop, legacy.state.inputs.crop)
    test_fixture_arrays_equal(adapted.state.inputs.soil, legacy.state.inputs.soil)
    test_fixture_arrays_equal(adapted.state.inputs.management, legacy.state.inputs.management)
    test_fixture_arrays_equal(adapted.state.output, legacy.state.output)
    test_fixture_arrays_equal(adapted.diagnostics, legacy.diagnostics)
    @test adapted.daily_weather.temp == legacy.daily_weather.temp
    @test adapted.daily_weather.prec == legacy.daily_weather.prec
    @test adapted.daily_weather.swr == legacy.daily_weather.swr
    @test adapted.daily_weather.lwr == legacy.daily_weather.lwr
    @test adapted.daily_weather.annual_co2 == legacy.daily_weather.annual_co2
end
