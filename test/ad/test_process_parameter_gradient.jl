using Agrocosm
using Enzyme
using Test

isdefined(@__MODULE__, :weather_attribution_fixture) ||
    include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

# `enzyme_process_parameter_gradient` is the process-attribution counterpart of
# the weather gradient: d(yield)/d(theta) for the parameters that define each
# stress mechanism. It is the quantity the paper's attribution figure reports,
# so the load-bearing test is agreement with a central difference through the
# ORDINARY production replay - not self-consistency of the AD path.
#
# Two hazards make this more than a formality.
#
#   * The sink and terminal-heat parameters sit in the tail of a 672-byte
#     CFTParameters, past Enzyme's default type-analysis offset limit, where a
#     reverse derivative comes back as an EXACT ZERO with no error and a
#     correct primal. Every assertion below therefore also checks the gradient
#     is nonzero where the finite difference says it should be.
#   * `_replace_cft_parameters` must be called INSIDE the differentiated block.
#     Called outside, theta is not on the tape and the gradient is
#     structurally zero - the same symptom as the offset truncation, from a
#     different cause.

"""Yield from the ordinary production replay with one CFT field set to `value`."""
function _replay_yield(case, name::Symbol, value; settings...)
    cft = Agrocosm.CFTParameters{eltype(case.forcing), Int32}(;
        (f => (f === name ? eltype(case.forcing)(value) : getfield(case.cft, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
    replay = Agrocosm.weather_harvest_replay(
        case.forcing, case.state, cft, case.parameters, case.climate, case.days,
        case.harvest_day; settings...,
    )
    return replay.yield
end

"""Central difference of yield with respect to one CFT field."""
function _finite_difference(case, name::Symbol, value, step; settings...)
    high = _replay_yield(case, name, value + step; settings...)
    low = _replay_yield(case, name, value - step; settings...)
    return (high - low) / (2 * step)
end

@testset "Process-parameter gradient matches finite differences" begin
    T = Float64
    steps = 24
    config = DiurnalConfig(; steps, shape = Agrocosm.diurnal_shape_code(:sinusoid))
    # The standalone exposure pass with organ temperature and both reproductive
    # mechanisms: the configuration the roadmap now calls production. Thresholds
    # are lowered because the fixture runs at a constant 19/25 C.
    settings = (; diurnal_config = nothing, organ_temperature = true,
                reproductive_sink = true, terminal_heat = true,
                heat_exposure_config = config, irrigation = false,
                nitrogen_limit_vcmax = true, crop_resp_fix = true)
    # The whole season, not a short window: grain set acts at fphu 0.45-0.70 and
    # grain filling at 0.70-0.95, so a window that does not span both leaves one
    # of the two parameters with no developmental weight and a genuinely zero
    # derivative - indistinguishable from the offset-truncation failure.
    case = weather_attribution_fixture(1; T, window_days = :season,
        heat_exposure_config = config, organ_temperature = true,
        reproductive_sink = true, terminal_heat = true,
        sterility_rate = 0.004, sterility_temperature = 19.0,
        filling_rate = 0.004, filling_temperature = 17.0)

    names = (:sterility_rate, :filling_rate)
    theta = T[getfield(case.cft, name) for name in names]
    result = enzyme_process_parameter_gradient(
        theta, names, case.forcing, case.state, case.cft, case.parameters,
        case.climate, case.days, case.harvest_day; block_days = 20, settings...,
    )

    @test result.primal ≈ result.production_yield rtol = 1e-3
    @test result.reverse_primal ≈ result.primal rtol = 1e-6
    @test length(result.gradient) == 2
    @test all(isfinite, result.gradient)
    @test result.parameter_names === names

    for (index, name) in enumerate(names)
        value = theta[index]
        # 10% of the value: the limit here is the model's own non-smoothness,
        # not arithmetic precision, exactly as for the weather gradient.
        fd = _finite_difference(case, name, value, value * 0.1; settings...)
        ad = result.gradient[index]
        @info "process-parameter gradient" name ad fd
        # A gradient of exactly zero is the offset-truncation signature, so it
        # is asserted against rather than tolerated.
        @test ad != 0
        @test sign(ad) == sign(fd)
        @test ad ≈ fd rtol = 0.15
    end

    # Both mechanisms reduce yield, so both derivatives are negative. If one
    # comes back positive the mechanism is wired backwards.
    @test all(<(0), result.gradient)
end

@testset "The parameter gradient refuses ill-posed requests" begin
    T = Float64
    config = DiurnalConfig(; steps = 24,
                             shape = Agrocosm.diurnal_shape_code(:sinusoid))
    case = weather_attribution_fixture(1; T, window_days = :season,
        heat_exposure_config = config, organ_temperature = true,
        reproductive_sink = true, sterility_rate = 0.004,
        sterility_temperature = 19.0)
    settings = (; organ_temperature = true, reproductive_sink = true,
                heat_exposure_config = config)

    # Length mismatch, unknown field, and a repeated name - the last one matters
    # because two gradients would silently sum into one slot.
    @test_throws DimensionMismatch enzyme_process_parameter_gradient(
        T[0.004], (:sterility_rate, :filling_rate), case.forcing, case.state,
        case.cft, case.parameters, case.climate, case.days, case.harvest_day;
        settings...)
    @test_throws ArgumentError enzyme_process_parameter_gradient(
        T[0.004], (:not_a_field,), case.forcing, case.state, case.cft,
        case.parameters, case.climate, case.days, case.harvest_day; settings...)
    @test_throws ArgumentError enzyme_process_parameter_gradient(
        T[0.004, 0.004], (:sterility_rate, :sterility_rate), case.forcing,
        case.state, case.cft, case.parameters, case.climate, case.days,
        case.harvest_day; settings...)
    @test_throws ArgumentError enzyme_process_parameter_gradient(
        T[0.004], (:sterility_rate,), case.forcing, case.state, case.cft,
        case.parameters, case.climate, case.days, case.harvest_day;
        block_days = 0, settings...)
end
