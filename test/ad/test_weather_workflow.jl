using Test, TOML
using Enzyme
include(joinpath(@__DIR__, "..", "..", "scripts", "attribution", "common.jl"))
include(joinpath(@__DIR__, "..", "..", "scripts", "attribution", "historical_helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "scripts", "attribution", "case_helpers.jl"))
include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

@testset "Isolated attribution workflow contracts" begin
    mktempdir() do scratch
        root = joinpath(scratch, "paper1_global_cell", "output", "default_parameter_attribution_v1")
        mkpath(root)
        config = Dict{String, Any}(
            "study" => Dict("results_root" => root, "parameter_source" => "version_defaults", "expected_commit" => "fixture"),
            "paths" => Dict{String, Any}(),
            "run" => Dict("backend" => "cpu", "nitrogen_limit_vcmax" => true, "crop_resp_fix" => true),
            "management" => Dict("mode" => "fixed", "fixed_year" => 2015, "fertilizer" => "yes", "manure" => true, "with_tillage" => true),
            "climate" => Dict("wind_variable" => "windspeed"),
            "cfts" => Dict("cft_ids" => [1, 2, 3, 9], "water_systems" => ["rainfed", "irrigated"]),
        )
        @test check_study_config(config; check_repository = false) == realpath(root)
        @test require_inside(joinpath(root, "nested", "result.nc"), root) == joinpath(realpath(root), "nested", "result.nc")
        @test_throws ErrorException require_inside(joinpath(root, "..", "old_baseline"), root)
        symlink(scratch, joinpath(root, "escape"))
        @test_throws ErrorException require_inside(joinpath(root, "escape", "result.nc"), root)
        invalid = deepcopy(config)
        invalid["run"]["nitrogen_limit_vcmax"] = false
        @test_throws ErrorException check_study_config(invalid; check_repository = false)
        for value in (0, -1, true, 1.5)
            invalid = deepcopy(config)
            invalid["attribution"] = Dict("window_days" => value)
            @test_throws ErrorException check_study_config(invalid; check_repository = false)
        end
        invalid = deepcopy(config)
        invalid["paths"]["pool_allocation"] = "old_warmup.nc"
        @test_throws ErrorException check_study_config(invalid; check_repository = false)
        path = joinpath(root, "manifest.toml")
        write_study_toml(path, Dict("source" => "fixture"), root)
        @test TOML.parsefile(path)["source"] == "fixture"
        @test_throws ErrorException write_study_toml(path, Dict(), root)
        @test length(file_sha(path)) == 64
        record = completion_record(config, path; stage = "fixture", warmup_years = 600)
        @test record["warmup_years"] == 600
        @test record["simulation_config_sha256"] == simulation_config_sha(config)
        @test record["environment_project_sha256"] == file_sha(Base.active_project())
        @test record["julia_version"] == string(VERSION)
        case_config = deepcopy(config)
        case_config["cases"] = [Dict("name" => "selected_later")]
        @test simulation_config_sha(case_config) == simulation_config_sha(config)
        case_config["climate"]["wind_variable"] = "changed"
        @test simulation_config_sha(case_config) != simulation_config_sha(config)

        # Non-contiguous canonical cells catch accidental compact-index/grid-index mixing.
        grid = AgrocosmData.GridIndex(Float32[-100, -99, -98], Float32[40, 41],
            Int32[-1 30; 10 -1; -1 20], Int32[10, 20, 30], Int32[2, 3, 1], Int32[1, 2, 2])
        selection = select_cells(grid, [1, 3])
        variables = production_output_variables()
        values = Dict{Symbol, Any}(Symbol(v.group, :_, v.field) => Float32[1, 2] for v in variables)
        chunks = [OutputChunk(1, :annual, [365], selection.cell_ids, values)]
        output_path = joinpath(root, "annual.nc")
        write_compact_history(output_path, grid, selection, chunks, [1901])
        NCDataset(output_path, "r") do dataset
            @test dataset["longitude"][:] == Float32[-99, -100]
            @test dataset["latitude"][:] == Float32[40, 41]
            @test dataset["cell_id"][:] == [10, 30]
            @test dataset["crop_yield"][:, 1] == [1, 2]
            @test dataset["year"][:] == [1901]
        end
        @test_throws ErrorException write_compact_history(output_path, grid, selection, chunks, [1901])

        fixture = weather_attribution_fixture(3; T = Float32)
        case = merge(fixture, (; anchor_year = 1901, selection = (cell_ids = [30],)))
        factual = weather_harvest_replay(case.forcing, case.state, case.cft,
            case.parameters, case.climate, case.days, case.harvest_day)
        case_path = joinpath(root, "reference_2001.nc")
        write_weather_case_output(case_path, case, (gradient = zeros(Float32, size(case.forcing)),),
            factual, case.forcing, factual)
        NCDataset(case_path, "r") do dataset
            @test dataset.attrib["cell_id"] == 30
            @test sum(dataset["ad_active_window"][:]) == length(case.days)
            @test dataset["factual_gpp"][factual.daily.day] == factual.daily.gpp
            @test dataset["factual_storage_carbon"].attrib["units"] == "gC m-2"
            @test isnan(dataset["factual_gpp"][1])
            @test dataset["reference_wind"][:] == case.forcing[:, 1, 5]
            @test dataset.attrib["schema_version"] == 2
            @test dataset["factual_harvest_event"][case.harvest_day] == 1
            @test dataset["factual_nitrogen_uptake"].attrib["diagnostic_kind"] == "daily_flux"
            @test dataset["factual_vcmax"].attrib["diagnostic_kind"] == "state_or_auxiliary"
            @test dataset["factual_soil_nitrate"][factual.daily.day] == factual.daily.soil_nitrate
            @test dataset["day_of_year"][365] == 365
        end
        reference = case.forcing .+ 1
        days = (first(case.days) + 2):(first(case.days) + 4)
        controls = weather_window_controls(case.forcing, reference, days)
        outside = setdiff(axes(controls, 1), days)
        @test controls[outside, :, :] == case.forcing[outside, :, :]
        @test controls[days, :, :] == reference[days, :, :]
        @test controls !== case.forcing && controls !== reference
        @test_throws ArgumentError weather_window_controls(case.forcing, reference, 0:1)
        @test_throws ArgumentError weather_window_controls(case.forcing, reference, 2:1)
        @test_throws DimensionMismatch weather_window_controls(case.forcing, reference[1:3, :, :], 1:2)
        # One fixed weather window, including all five channels, must not affect
        # the causal prefix. No-op restoration must reproduce all diagnostics.
        noop = weather_harvest_replay(weather_window_controls(case.forcing, case.forcing, days),
            case.state, case.cft, case.parameters, case.climate, case.days, case.harvest_day)
        @test noop.daily == factual.daily && noop.yield == factual.yield
        replay = weather_harvest_replay(controls, case.state, case.cft, case.parameters,
            case.climate, case.days, case.harvest_day; replay_end_day = size(controls, 1))
        for name in keys(factual.daily)
            @test getproperty(replay.daily, name)[1:2] == getproperty(factual.daily, name)[1:2]
        end
        window_path = joinpath(root, "window.nc")
        record = write_weather_counterfactual(window_path, case,
            (gradient = zeros(Float32, size(case.forcing)),), factual, controls, replay;
            replacement_days = days, counterfactual_kind = "window_restore", reference_year = 1902)
        @test record["yield_change"] == replay.yield - factual.yield
        @test record["output_sha256"] == file_sha(window_path)
        NCDataset(window_path, "r") do dataset
            @test findall(==(1), dataset["weather_replacement_window"][:]) == collect(days)
            @test dataset.attrib["reference_year"] == 1902
            @test dataset.attrib["counterfactual_kind"] == "window_restore"
        end
        @test_throws ErrorException write_weather_case_output(window_path, case,
            (gradient = zeros(Float32, size(case.forcing)),), factual, controls, replay)
    end
end
