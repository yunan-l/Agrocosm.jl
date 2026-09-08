using Agrocosm, Enzyme, Test
include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

function weather_kernel_probe(controls, weather)
    Agrocosm.apply_weather_forcing!(weather, controls, 2)
    return sum(weather.temp) + 2sum(weather.prec) + 3sum(weather.swr) +
        4sum(weather.lwr) + 5sum(weather.wind)
end

@testset "All five weather channels remain active" begin
    controls = ones(3, 2, 5)
    gradient = zeros(size(controls))
    weather = Agrocosm.init_weather(Float64, 2, identity)
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal), weather_kernel_probe,
        Enzyme.Duplicated(controls, gradient), Enzyme.Duplicated(weather, Enzyme.make_zero(weather)))
    for variable in 1:5
        @test gradient[2, :, variable] == fill(variable, 2)
    end
    @test all(iszero, gradient[[1, 3], :, :])
end

@testset "Weather replay and event boundaries" begin
    for T in (Float32, Float64), id in (1, 3)
        case = weather_attribution_fixture(id; T)
        reference = weather_harvest_replay(case.forcing, case.state, case.cft,
            case.parameters, case.climate, case.days, case.harvest_day)
        @test reference.schedule_matches && !reference.failed
        @test length(reference.daily.day) == length(case.days) + 1
        @test all(isfinite, reference.daily.gpp)
        @test all(isfinite, reference.daily.et)
        @test reference.daily.harvest_event == [zeros(Int8, length(case.days)); 1]
        @test all(name -> all(isfinite, getproperty(reference.daily, name)), keys(reference.daily))
        growth = 1:length(case.days)
        @test reference.daily.npp[growth] ≈ reference.daily.gpp[growth] .-
            reference.daily.respiration[growth] .- reference.daily.biological_fixation_cost[growth]
        @test reference.daily.biomass_carbon[growth] ≈ reference.daily.leaf_carbon[growth] .+
            reference.daily.root_carbon[growth] .+ reference.daily.mobile_carbon[growth] .+
            reference.daily.storage_carbon[growth]
        @test all(x -> 0 <= x <= 1, reference.daily.water_sufficiency[growth])
        @test all(x -> 0 <= x <= 1, reference.daily.nitrogen_limitation[growth])
        @test weather_forcing(case.climate) == case.forcing
    end
    # Exactly 600 PHU in constant-temperature maize puts the factual trajectory
    # on an event threshold. Its temperature finite difference is not a smooth
    # derivative of a fixed-harvest objective.
    case = weather_attribution_fixture(3; phu = 600)
    controls = copy(case.forcing)
    controls[first(case.days):(first(case.days) + 2), 1, 1] .-= 0.01
    delayed = weather_harvest_replay(controls, case.state, case.cft, case.parameters,
        case.climate, case.days, case.harvest_day)
    @test !delayed.schedule_matches
    @test isempty(delayed.harvest_days)
    @test isnan(delayed.yield) # Unharvested is not silently reported as zero.
    extended = weather_harvest_replay(controls, case.state, case.cft, case.parameters,
        case.climate, case.days, case.harvest_day; replay_end_day = case.harvest_day + 3)
    @test extended.harvest_days == [case.harvest_day + 1]
    @test isfinite(extended.yield) && extended.yield > 0
    @test !extended.schedule_matches
    @test_throws ArgumentError enzyme_weather_harvest_gradient(controls, case.state,
        case.cft, case.parameters, case.climate, case.days, case.harvest_day)
    weather = Agrocosm.init_weather(Float64, 1, identity)
    @test_throws BoundsError Agrocosm.apply_weather_forcing!(weather, case.forcing, 0)
    @test_throws DimensionMismatch Agrocosm.apply_weather_forcing!(weather, zeros(365, 1, 4), 1)
end
