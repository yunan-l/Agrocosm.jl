#!/usr/bin/env julia
include(joinpath(@__DIR__, "common.jl"))
length(ARGS) == 2 || error("usage: run_attribution.jl STUDY.toml CASE_NAME")
config_path = abspath(ARGS[1])
config = TOML.parsefile(config_path)
root = check_study_config(config)
specification = only(filter(item -> item["name"] == ARGS[2], get(config, "cases", [])))
check_patch(config, specification["cft_id"], specification["water_system"])
haskey(ENV, "SLURM_JOB_ID") || error("historical case preparation is server/Slurm-only")

using Enzyme
using JLD2
include(joinpath(@__DIR__, "historical_helpers.jl"))
include(joinpath(@__DIR__, "case_helpers.jl"))
check_loaded_attribution_model()
directory = claim_stage_directory(config, "attribution", specification["name"])
case = prepare_historical_case(config, config_path, specification)
irrigation = specification["water_system"] == "irrigated"
jldsave(joinpath(directory, "prepared_case.jld2");
    initial_state = case.state, cft = case.cft, parameters = case.parameters,
    climate = case.climate, days = case.days, harvest_day = case.harvest_day,
    specification, source_commit = config["study"]["expected_commit"])
result = enzyme_weather_harvest_gradient(case.forcing, case.state, case.cft,
    case.parameters, case.climate, case.days, case.harvest_day;
    irrigation, block_days = config["attribution"]["block_days"])
factual = weather_harvest_replay(case.forcing, case.state, case.cft,
    case.parameters, case.climate, case.days, case.harvest_day; irrigation)
T = eltype(case.forcing)
direction = zeros(T, size(case.forcing))
direction[first(case.days):min(first(case.days) + 2, last(case.days)), 1, 1] .= one(T)
projection = sum(result.gradient .* direction)
validation = Dict{String, Any}("status" => "event_shift_at_both_steps",
    "probe" => "temperature +1 C in first three post-sowing days", "ad_projection" => projection)
for step in T[0.05, 0.01]
    plus = weather_harvest_replay(case.forcing .+ step .* direction, case.state, case.cft,
        case.parameters, case.climate, case.days, case.harvest_day; irrigation)
    minus = weather_harvest_replay(case.forcing .- step .* direction, case.state, case.cft,
        case.parameters, case.climate, case.days, case.harvest_day; irrigation)
    if plus.schedule_matches && minus.schedule_matches && !plus.failed && !minus.failed
        finite_difference = (plus.yield - minus.yield) / (2step)
        isapprox(projection, finite_difference; rtol = 0.03, atol = 2e-5) ||
            error("case weather gradient fails ordinary-production finite difference: $projection vs $finite_difference")
        merge!(validation, Dict("status" => "passed", "step" => step,
            "finite_difference" => finite_difference))
        break
    end
end
validation["status"] == "passed" || @warn "Finite difference crosses a harvest/failure boundary; gradient is conditional only" validation
references = specification["reference_years"]
!isempty(references) && allunique(references) || error("choose unique, explicit reference weather years")
specification["event_year"] in references && error("reference years must differ from the event year")
records = Dict{String, Any}[]
window_days = get(config["attribution"], "window_days", 14)
for reference_year in references
    controls = reference_weather(config, case, reference_year)
    counterfactual = weather_harvest_replay(controls, case.state, case.cft,
        case.parameters, case.climate, case.days, case.harvest_day;
        irrigation, replay_end_day = size(controls, 1))
    path = joinpath(directory, "reference_$(reference_year).nc")
    record = write_weather_counterfactual(path, case, result, factual, controls, counterfactual;
        replacement_days = first(case.days):size(controls, 1),
        counterfactual_kind = "full_reference", reference_year)
    windows = Dict{String, Any}[]
    # Disjoint, predeclared calendar windows; not selected from the largest
    # observed response. Reuse the factual gradient, then verify each finite
    # change with production. Do not sum these nonlinear single-window effects.
    for start_day in first(case.days):window_days:last(case.days)
        days = start_day:min(start_day + window_days - 1, last(case.days))
        window_controls = weather_window_controls(case.forcing, controls, days)
        replay = weather_harvest_replay(window_controls, case.state, case.cft,
            case.parameters, case.climate, case.days, case.harvest_day;
            irrigation, replay_end_day = size(controls, 1))
        window_path = joinpath(directory, "reference_$(reference_year)_window_$(first(days))_$(last(days)).nc")
        push!(windows, write_weather_counterfactual(window_path, case, result, factual, window_controls, replay;
            replacement_days = days, counterfactual_kind = "window_restore", reference_year))
    end
    record["windows"] = windows
    push!(records, record)
end
write_study_toml(joinpath(directory, "study_complete.toml"),
    completion_record(config, config_path; stage = "attribution", case = specification,
        production_yield = factual.yield, historical_yield = case.historical_yield,
        ad_primal = result.primal, reverse_primal = result.reverse_primal,
        directional_validation = validation,
        harvest_day = case.harvest_day, anchor_year = case.anchor_year,
        attribution_scope = "conditional_on_post_sowing_state_and_fixed_harvest",
        window_days, process_scope = "paired weather-window process responses; not additive process attribution",
        allocation_sha256 = file_sha(case.allocation_file), references = records), root)
@info "Attribution case complete" directory yield = factual.yield
