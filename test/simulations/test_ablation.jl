using Agrocosm
using Test

# The ablation ladder is a scientific object: the paper reports it as the
# sequence of structural deficiencies it isolates, so the registry and the code
# must not be able to drift apart. These tests are the mechanism that stops
# them. The interesting ones are not "does the table say what I typed" but
# "does the constructor actually enforce what the table claims" and "does every
# rung the table offers actually build".

_configuration(; kwargs...) = Agrocosm.SimulationConfiguration(
    Float32, identity, 10, [1], [1]; kwargs...,
)

@testset "Ablation ladder is consistent with the configuration contract" begin
    @test first(ablation_rungs()) === :daily
    @test ablation_rungs() == (:daily, :subdaily, :organ_temperature, :reproductive_sink)
    @test length(ABLATION_LADDER) == length(ablation_rungs()) - 1

    # Every step must name a real field, or `ablation_configuration` would
    # silently produce a keyword the constructor ignores.
    fields = fieldnames(Agrocosm.SimulationConfiguration)
    for step in ABLATION_LADDER
        @test step.field in fields
        @test !isempty(step.summary)
        @test !isempty(step.retreat)
    end

    # Every declared prerequisite must name a step strictly below the one that
    # declares it, so the cumulative sequence is well founded. Note this is
    # weaker than "the previous step": the sink's prerequisite is the sub-daily
    # loop, not organ temperature, which is exactly what makes the air-driven
    # sink cell expressible.
    @test ABLATION_LADDER[1].requires === nothing
    names = map(step -> step.name, ABLATION_LADDER)
    for index in 2:length(ABLATION_LADDER)
        required = ABLATION_LADDER[index].requires
        @test required !== nothing
        @test findfirst(==(required), names) < index
    end
    @test ablation_step(:reproductive_sink).requires === :subdaily

    @test_throws ArgumentError ablation_step(:nonexistent)
    @test ablation_step(:subdaily).field === :subdaily_photosynthesis
end

@testset "Ablation prerequisites are enforced, not just documented" begin
    # This is the drift guard. If someone relaxes the validation in
    # `SimulationConfiguration`, a registry entry claiming a prerequisite becomes
    # untrue and this fails - even though the registry still reads fine.
    for index in 2:length(ABLATION_LADDER)
        step = ABLATION_LADDER[index]
        prerequisite = ablation_step(step.requires)
        @test_throws ArgumentError _configuration(;
            step.field => true, prerequisite.field => false,
        )
    end
end

@testset "The air-driven sink cell is expressible and off the ladder" begin
    # The comparison that separates the sink mechanism from the leaf-air
    # departure triggering it: same sink, same threshold, same sub-daily
    # resolution, air temperature instead of leaf temperature.
    settings = ablation_air_driven_sink_configuration(; subdaily_steps = 24)
    @test settings.subdaily_photosynthesis
    @test !settings.organ_temperature
    @test settings.reproductive_sink
    configuration = _configuration(; settings...)
    @test configuration.reproductive_sink
    @test !configuration.organ_temperature

    # It is deliberately NOT a rung, so no rung reproduces it. Every ladder rung
    # with the sink on also has organ temperature on.
    for rung in ablation_rungs()
        rung_settings = ablation_configuration(rung)
        rung_settings.reproductive_sink && @test rung_settings.organ_temperature
        @test rung_settings != settings
    end

    # And the sink still cannot stand without the sub-daily loop, which is the
    # prerequisite that did not get relaxed.
    @test_throws ArgumentError _configuration(;
        reproductive_sink = true, subdaily_photosynthesis = false,
    )
    for field in (:subdaily_photosynthesis, :organ_temperature, :reproductive_sink)
        @test_throws ArgumentError ablation_air_driven_sink_configuration(;
            field => true,
        )
    end
end

@testset "Each ablation rung is a complete, buildable structure statement" begin
    owned = map(step -> step.field, ABLATION_LADDER)
    for (position, rung) in enumerate(ablation_rungs())
        settings = ablation_configuration(rung)
        # Complete: every switch the ladder owns is stated, so a rung is not a
        # patch on whatever the defaults happen to be. Fig 1's "identical
        # parameters" control depends on this.
        @test Set(keys(settings)) == Set(owned)
        # On up to the rung, off above it. `:daily` is position 1, below every
        # step, so nothing is on there.
        for (index, step) in enumerate(ABLATION_LADDER)
            @test settings[step.field] == (index <= position - 1)
        end
        # And it has to actually build.
        configuration = _configuration(; settings...)
        for step in ABLATION_LADDER
            @test getfield(configuration, step.field) == settings[step.field]
        end
    end

    @test all(iszero, values(ablation_configuration(:daily)))
    @test all(values(ablation_configuration(last(ablation_rungs()))))
    @test_throws ArgumentError ablation_configuration(:hourly)
end

@testset "Ablation rungs carry shared settings and refuse to be overridden" begin
    settings = ablation_configuration(:subdaily; subdaily_steps = 48,
                                      diurnal_shape = :daytime_neutral,
                                      irrigation = true)
    @test settings.subdaily_photosynthesis
    @test settings.subdaily_steps == 48
    @test settings.diurnal_shape === :daytime_neutral
    @test settings.irrigation
    configuration = _configuration(; settings...)
    @test configuration.subdaily_steps == 48
    @test configuration.diurnal_shape === :daytime_neutral

    # Overriding a switch the ladder owns would move the run off the ladder
    # while still labelling it with a rung name, which is the one mistake this
    # interface exists to prevent.
    for step in ABLATION_LADDER
        @test_throws ArgumentError ablation_configuration(:daily; step.field => true)
        @test_throws ArgumentError ablation_configuration(:subdaily; step.field => false)
    end

    # The whole ladder in one call, sharing settings, which is what a driver
    # needs to keep the comparison controlled.
    ladder = ablation_ladder_settings(; subdaily_steps = 24)
    @test map(first, ladder) == ablation_rungs()
    @test all(pair -> last(pair).subdaily_steps == 24, ladder)
    @test !last(ladder[1]).subdaily_photosynthesis
    @test last(ladder[end]).reproductive_sink
end
