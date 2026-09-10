using Agrocosm
using Test

include(joinpath(@__DIR__, "..", "..", "scripts", "run_global_wheat_cpu.jl"))

# THE test that was missing, and whose absence let a gap survive that would have
# consumed a server allocation.
#
# The global entry point forwarded `crop_resp_fix`, `nitrogen_limit_vcmax`,
# `sowing_mode` and `freeze_vernalization_requirement` to
# `initialize_simulation` - and nothing else. None of the ten process switches
# this project added reached the constructor, so every global run executed rung
# zero regardless of what its config said, completed normally, and wrote
# plausible output.
#
# The failure was silent BY DESIGN and correctly so: the model refuses an
# inconsistent configuration, but a configuration with everything off is
# perfectly valid - it is the paper's control. Nothing could have detected this
# except a test that asserts a named configuration actually switches something
# on, which is what this file is.

@testset "Every named configuration resolves to distinct switches" begin
    names = (:daily, :subdaily, :organ_temperature, :reproductive_sink,
             :tmax_sink, :daily_sink, :daily_sink_air, :terminal_only, :production)
    resolved = Dict(name => process_settings(
        Dict("processes" => Dict("configuration" => String(name)))) for name in names)

    # Distinctness is the property that matters: if two names resolved to the
    # same switches, the ablation table would have two rows measuring one thing.
    switch_view(settings) = (
        get(settings, :subdaily_photosynthesis, false),
        get(settings, :subdaily_heat_exposure, false),
        get(settings, :daily_statistic_exposure, false),
        get(settings, :organ_temperature, false),
        get(settings, :reproductive_sink, false),
        get(settings, :terminal_heat, false),
        get(settings, :water_sterility, false),
        get(settings, :water_filling, false),
    )
    views = Dict(name => switch_view(settings) for (name, settings) in resolved)
    @test length(unique(values(views))) == length(names)

    # Rung zero must have nothing on, and be what an absent section gives.
    @test all(iszero, views[:daily])
    @test switch_view(process_settings(Dict{String, Any}())) == views[:daily]

    # The production configuration must have every mechanism on. This is the
    # assertion the gap would have failed.
    @test views[:production] == (false, true, false, true, true, true, true, true)

    # And every named configuration must build.
    for (name, settings) in resolved
        configuration = Agrocosm.SimulationConfiguration(
            Float32, identity, 10, [1], [1]; settings...,
        )
        @test configuration isa Agrocosm.SimulationConfiguration
        # Exactly one writer for the exposure fields, or the field carries
        # whichever kernel ran last.
        writers = count((configuration.subdaily_photosynthesis,
                         configuration.subdaily_heat_exposure,
                         configuration.daily_statistic_exposure))
        @test writers <= 1
        # A heat mechanism with no writer is inert, which is never intended.
        (configuration.reproductive_sink || configuration.terminal_heat) &&
            @test writers == 1
    end

    @test_throws ArgumentError process_settings(
        Dict("processes" => Dict("configuration" => "not_a_configuration")))
end

@testset "The rate ray scales all four reproductive rates together" begin
    # The bound was taken by scaling all four at once, so a config that scales a
    # subset is not on the ray the bound lives on.
    base = Agrocosm.cft1
    for scale in (0.0, 0.5, 2.0)
        scaled = scaled_cft(base, scale)
        for field in (:sterility_rate, :filling_rate,
                      :water_sterility_rate, :water_filling_rate)
            @test getfield(scaled, field) ≈ scale * getfield(base, field)
        end
        # Nothing else may move: a rate sweep that also changed a threshold
        # would confound the two.
        for field in (:sterility_temperature, :filling_temperature, :hiopt,
                      :flowering_start, :flowering_end, :laimax)
            @test getfield(scaled, field) == getfield(base, field)
        end
    end
    # Scale one is the identity, exactly, so the default path is untouched.
    @test scaled_cft(base, 1) === base
end

@testset "Shared sub-daily settings reach the configuration" begin
    settings = process_settings(Dict("processes" => Dict(
        "configuration" => "production", "subdaily_steps" => 48,
        "diurnal_shape" => "daytime_neutral")))
    @test settings.subdaily_steps == 48
    @test settings.diurnal_shape === :daytime_neutral
    configuration = Agrocosm.SimulationConfiguration(
        Float32, identity, 10, [1], [1]; settings...,
    )
    @test configuration.subdaily_steps == 48
    @test Agrocosm.heat_exposure_configuration(configuration) !== nothing
end
