using Agrocosm
using Enzyme
using Test

isdefined(@__MODULE__, :weather_attribution_fixture) ||
    include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

# d(yield)/d(heat_day_rate) and d(yield)/d(cold_night_rate) through the
# differentiated path.
#
# This file exists because both were structurally ZERO and nothing noticed.
# `_enzyme_continuous_transition!` inlines the daily sequence rather than calling
# `_daily_crop!`, and it carried `reproductive_sink` and `terminal_heat` but not
# `anthesis_heat`, `water_sterility`, `water_filling` or `cold_sterility`. A
# mechanism absent from that function contributes nothing to the gradient the
# paper's attribution figure reports, with a correct primal and no error - the
# same silent-zero signature `test_process_parameter_gradient.jl` warns about
# for the type-analysis offset limit, from a different cause.
#
# `cold_night_rate` is now the LAST field of `CFTParameters`, deepest into the
# tail where that offset truncation bites, so the nonzero assertions here carry
# more weight than usual.
#
# `water_sterility` and `water_filling` are STILL missing from that function.
# They are off-ladder and older, and wiring them is a separate change; this
# comment is the record that it is outstanding.

"""Yield from the ordinary production replay with one CFT field set to `value`."""
function _threshold_replay_yield(case, name::Symbol, value; settings...)
    cft = Agrocosm.CFTParameters{eltype(case.forcing), Int32}(;
        (f => (f === name ? eltype(case.forcing)(value) : getfield(case.cft, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
    replay = Agrocosm.weather_harvest_replay(
        case.forcing, case.state, cft, case.parameters, case.climate, case.days,
        case.harvest_day; settings...,
    )
    return replay.yield
end

function _threshold_finite_difference(case, name::Symbol, value, step; settings...)
    high = _threshold_replay_yield(case, name, value + step; settings...)
    low = _threshold_replay_yield(case, name, value - step; settings...)
    return (high - low) / (2 * step)
end

@testset "Absolute-threshold rates reach the gradient" begin
    T = Float64
    # The fixture holds 19/25 C with a 10 C range, so the daily maximum is 24 C
    # and the daily minimum 14 C. Thresholds are placed inside that band for the
    # same reason the sink's are: the real 38 C and 17 C are never reached here.
    settings = (; anthesis_heat = true, cold_sterility = true,
                irrigation = false, nitrogen_limit_vcmax = true,
                crop_resp_fix = true)
    case = weather_attribution_fixture(1; T, window_days = :season,
        anthesis_heat = true, cold_sterility = true,
        heat_day_temperature = 22.0, heat_day_rate = 0.004,
        cold_night_temperature = 17.0, cold_night_rate = 0.004)

    # Both must actually bite in this fixture, or every assertion below is
    # vacuous: a mechanism that never fires has a genuinely zero derivative.
    baseline = _threshold_replay_yield(case, :heat_day_rate, 0.0; settings...)
    @test _threshold_replay_yield(case, :heat_day_rate, 0.01; settings...) < baseline
    cold_baseline = _threshold_replay_yield(case, :cold_night_rate, 0.0; settings...)
    @test _threshold_replay_yield(case, :cold_night_rate, 0.01; settings...) <
          cold_baseline

    names = (:heat_day_rate, :cold_night_rate)
    theta = T[getfield(case.cft, name) for name in names]
    result = enzyme_process_parameter_gradient(
        theta, names, case.forcing, case.state, case.cft, case.parameters,
        case.climate, case.days, case.harvest_day; block_days = 20, settings...,
    )

    @test result.primal ≈ result.production_yield rtol = 1e-3
    @test length(result.gradient) == 2
    @test all(isfinite, result.gradient)

    for (index, name) in enumerate(names)
        value = theta[index]
        fd = _threshold_finite_difference(case, name, value, value * 0.1; settings...)
        ad = result.gradient[index]
        @info "absolute-threshold gradient" name ad fd
        # Exactly zero is the signature this file was written to catch.
        @test ad != 0
        @test sign(ad) == sign(fd)
        @test ad ≈ fd rtol = 0.15
    end
    # Both mechanisms destroy grain set, so both derivatives are negative.
    @test all(<(0), result.gradient)
end

@testset "The positional order of the AD transition flags is preserved" begin
    # `_enzyme_continuous_transition!` takes its mechanism flags POSITIONALLY
    # with defaults, and `weather_attribution.jl` passes the whole list
    # positionally. Adding the two new flags anywhere but the END silently
    # rebinds `terminal_heat` to `anthesis_heat` at that one call site, with no
    # error and a plausible result. Assert the order rather than the convention.
    extension = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)
    @test extension !== nothing
    # Optional positional arguments expand into one method per arity - fifteen
    # of them here - so take the widest rather than the first, which is an
    # eight-argument stub.
    candidates = collect(methods(extension._enzyme_continuous_transition!))
    method = argmax(m -> m.nargs, candidates)
    names = Base.method_argnames(method)
    @test names[end] === :cold_sterility
    @test names[end - 1] === :anthesis_heat
    @test names[end - 2] === :terminal_heat
end
