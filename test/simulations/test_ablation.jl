using Agrocosm
using Test

# The ablation ladder is a scientific object: the paper reports it as the
# sequence of structural deficiencies it isolates, so the registry and the code
# must not be able to drift apart. These tests are the mechanism that stops
# them. The interesting ones are not "does the table say what I typed" but
# "does the constructor actually enforce what the table claims" and "does every
# rung the table offers actually build".


"""Drive `ablation_metrics`' reduction on synthetic annual output.

Mirrors the real reduction rather than re-deriving it, so the test fails if the
reduction changes. Built as a minimal stand-in for a finished `CropSimulation`'s
annual rows: `yield` and `harvest_aboveground_carbon` in gC, and the annual
`harvest_date`, which is nonzero for any season that ended - destroyed or not.
"""
function _fake_metrics(season_yield, season_above, harvest_date)
    scale = Agrocosm.GRAMS_PER_M2_TO_TONNES_PER_HECTARE /
            Agrocosm.CARBON_FRACTION_OF_DRY_MATTER
    completed = count(>(0), harvest_date)
    with_grain = count(>(0), season_yield)
    ended = [(g, a) for (g, a) in zip(season_yield, season_above) if a > 0]
    total = sum(filter(>(0), season_yield); init = 0.0) * scale
    return (
        seasons = completed,
        seasons_with_grain = with_grain,
        seasons_destroyed = completed - with_grain,
        yield_total = total,
        yield_per_season = completed == 0 ? 0.0 : total / completed,
        harvest_index = isempty(ended) ? 0.0 :
            sum(g / a for (g, a) in ended) / length(ended),
    )
end

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

@testset "The daily-assimilation sink cell is expressible and off the ladder" begin
    # The configuration today's separability measurement argues for: the sink
    # and leaf temperature without the sub-daily assimilation loop, so the
    # calibrated daily kernel grows the canopy that the exposure integral then
    # reads. See docs/07_ablation_framework.md.
    settings = ablation_daily_assimilation_sink_configuration()
    @test !settings.subdaily_photosynthesis
    @test settings.subdaily_heat_exposure
    @test settings.organ_temperature
    @test settings.reproductive_sink
    configuration = _configuration(; settings...)
    @test !configuration.subdaily_photosynthesis
    @test configuration.subdaily_heat_exposure
    @test configuration.reproductive_sink
    # Exactly one loop is live, which is what makes the exposure field have one
    # writer.
    @test Agrocosm.diurnal_configuration(configuration) === nothing
    @test Agrocosm.heat_exposure_configuration(configuration) !== nothing

    # The leaf/air comparison has to survive into this architecture too, since
    # air temperature loses most of the exposure hours.
    air = ablation_daily_assimilation_sink_configuration(; organ_temperature = false)
    @test !air.organ_temperature
    @test air.subdaily_heat_exposure
    @test _configuration(; air...).reproductive_sink

    # Not a rung: no point on the ladder reproduces it, because every rung with
    # the sink on also has the sub-daily assimilation loop on.
    for rung in ablation_rungs()
        rung_settings = ablation_configuration(rung)
        @test !haskey(rung_settings, :subdaily_heat_exposure)
        rung_settings.reproductive_sink && @test rung_settings.subdaily_photosynthesis
    end

    # A fixed comparison cell apart from `organ_temperature`, so the switches it
    # owns cannot be overridden into something still carrying its name.
    for field in (:subdaily_photosynthesis, :subdaily_heat_exposure, :reproductive_sink)
        @test_throws ArgumentError ablation_daily_assimilation_sink_configuration(;
            field => true,
        )
    end
end

@testset "The daily-statistic floor cell is expressible and off the ladder" begin
    # The cell that keeps "a daily model responds exactly zero" from being a
    # definition: tasmax rises when the range widens, so a criterion built on it
    # is not inert, and the cell measures what it recovers.
    settings = ablation_daily_statistic_sink_configuration()
    @test settings.daily_statistic_exposure
    @test !settings.subdaily_photosynthesis
    @test !settings.subdaily_heat_exposure
    @test !settings.organ_temperature
    @test settings.reproductive_sink
    configuration = _configuration(; settings...)
    @test Agrocosm.daily_statistic_exposure_enabled(configuration)
    @test Agrocosm.diurnal_configuration(configuration) === nothing
    @test Agrocosm.heat_exposure_configuration(configuration) === nothing

    # No rung reproduces it, and it owns every switch it sets.
    for rung in ablation_rungs()
        @test !haskey(ablation_configuration(rung), :daily_statistic_exposure)
    end
    for field in (:subdaily_photosynthesis, :subdaily_heat_exposure,
                  :daily_statistic_exposure, :organ_temperature, :reproductive_sink)
        @test_throws ArgumentError ablation_daily_statistic_sink_configuration(;
            field => true,
        )
    end

    # The closed form has no sub-steps, so it cannot host a canopy energy
    # balance; asking for one is a contradiction rather than a silent no-op.
    @test_throws ArgumentError _configuration(;
        daily_statistic_exposure = true, organ_temperature = true,
    )
end

@testset "heat_exposure_hours has exactly one writer" begin
    # Two writers would leave the field carrying whichever kernel ran last, so
    # every ablation cell reading it would silently stop measuring what it
    # claims. The constructor is where that has to fail.
    for (a, b) in ((:subdaily_photosynthesis, :subdaily_heat_exposure),
                   (:subdaily_photosynthesis, :daily_statistic_exposure),
                   (:subdaily_heat_exposure, :daily_statistic_exposure))
        @test_throws ArgumentError _configuration(; a => true, b => true)
    end
    @test_throws ArgumentError _configuration(;
        subdaily_photosynthesis = true, subdaily_heat_exposure = true,
        daily_statistic_exposure = true,
    )
    # All three satisfy the sink; only the two sub-daily loops satisfy organ
    # temperature, and nothing is satisfied by no source at all.
    for host in (:subdaily_photosynthesis, :subdaily_heat_exposure,
                 :daily_statistic_exposure)
        @test _configuration(; host => true, reproductive_sink = true).reproductive_sink
    end
    for host in (:subdaily_photosynthesis, :subdaily_heat_exposure)
        @test _configuration(; host => true, reproductive_sink = true).reproductive_sink
        @test _configuration(; host => true, organ_temperature = true).organ_temperature
    end
    @test_throws ArgumentError _configuration(; reproductive_sink = true)
    @test_throws ArgumentError _configuration(; organ_temperature = true)
    # The capacity solve belongs to the assimilation light course, so the
    # standalone pass does not host it.
    @test_throws ArgumentError _configuration(;
        subdaily_heat_exposure = true, subdaily_capacity_optimum = true,
    )
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

@testset "A destroyed season lowers the ablation metric" begin
    # The defect this guards against is subtle and was live: dividing yield by
    # the count of seasons that YIELDED removes a destroyed season from the
    # numerator and the denominator together, so an event that annihilates a
    # crop reports no change at all. On an extreme-event ladder that is not a
    # rough edge, it is the metric being blind to the outcome it exists to
    # measure. Two seasons, one destroyed, must read lower than two intact
    # ones.
    intact = _fake_metrics([200.0, 200.0], [400.0, 400.0], Int32[150, 515])
    ruined = _fake_metrics([200.0, 0.0], [400.0, 400.0], Int32[150, 515])
    absent = _fake_metrics([200.0, 0.0], [400.0, 0.0], Int32[150, 0])

    @test intact.seasons == 2 && intact.seasons_with_grain == 2
    @test intact.seasons_destroyed == 0
    @test ruined.seasons == 2 && ruined.seasons_with_grain == 1
    @test ruined.seasons_destroyed == 1
    @test ruined.yield_per_season < intact.yield_per_season
    @test ruined.yield_per_season ≈ intact.yield_per_season / 2
    # And the harvest index registers it too: a season with above-ground carbon
    # and no grain contributes a genuine zero.
    @test ruined.harvest_index ≈ intact.harvest_index / 2

    # A row where no season happened at all is not a destroyed season: no
    # harvest date, no above-ground carbon, so it must not dilute either mean.
    @test absent.seasons == 1
    @test absent.seasons_destroyed == 0
    @test absent.yield_per_season ≈ intact.yield_per_season
    @test absent.harvest_index ≈ intact.harvest_index
end
