using Test
using Agrocosm

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Temperature stress over the day's radiation-weighted temperature course
# instead of at the daily mean.
#
# This is a change to the ARGUMENT of a response function the model already has,
# not a new damage term, and the two contracts that matter are that it is bitwise
# the old form when off, and that when on it reproduces the numbers
# `tools/quantify_aggregation_bias.py` reports - that tool is what measured the
# effect before the kernel existed, so agreement is a cross-check between an
# independent implementation and this one.

const TD = Float64

function stress_state(temp, range; steps, cft = Agrocosm.cft1, daylength = 14.0)
    crop = test_model_state(Agrocosm.init_crop(TD, 1, identity))
    crop.auxiliary.pet.daylength .= TD(daylength)
    Agrocosm.temp_stress(
        Agrocosm.convert_precision(TD, cft), crop.auxiliary.pet, crop,
        fill(TD(temp), 1);
        diurnal_range = steps > 0 ? fill(TD(range), 1) : nothing,
        diurnal_steps = steps,
        diurnal_shape = Agrocosm.diurnal_shape_code(:sinusoid),
    )
    return Agrocosm.crop_photosynthesis_auxiliary(crop).temperature_stress[1]
end

@testset "steps zero is bitwise the daily-mean form" begin
    for temp in TD[10.0, 20.0, 28.0, 34.0, 40.0], range in TD[0.0, 8.0, 16.0]
        @test stress_state(temp, range; steps = 0) === stress_state(temp, 0.0; steps = 0)
    end
    # And a range with no steps is ignored rather than half-applied.
    @test stress_state(34.0, 12.0; steps = 0) === stress_state(34.0, 0.0; steps = 0)
end

@testset "a flat day reproduces the daily mean" begin
    # Zero amplitude means every sub-step is at the mean, and the radiation
    # weights sum to one, so the weighted mean of a constant is that constant.
    for temp in TD[15.0, 25.0, 34.0]
        @test stress_state(temp, 0.0; steps = 24) ≈ stress_state(temp, 0.0; steps = 0) rtol = 1e-12
    end
end

@testset "it reproduces the independently computed aggregation bias" begin
    # `quantify_aggregation_bias.py`, sinusoid, daylength 14 h, DTR 12 C.
    # Its numbers were measured before this kernel was written.
    for (cft, mean_temp, expected) in (
        (Agrocosm.cft1, 34.0, 0.3079), (Agrocosm.cft1, 36.0, 0.1521),
        (Agrocosm.cft3, 34.0, 0.6515), (Agrocosm.cft3, 36.0, 0.3810),
        (Agrocosm.cft2, 34.0, 0.9995),
    )
        @test stress_state(mean_temp, 12.0; steps = 48, cft) ≈ expected atol = 5e-3
    end
    # Rice is the crop the change cannot touch: its photosynthesis optimum runs
    # to 45 C, so the course never leaves the plateau.
    @test stress_state(34.0, 12.0; steps = 48, cft = Agrocosm.cft2) >
          stress_state(34.0, 12.0; steps = 48, cft = Agrocosm.cft1)
end

@testset "the course costs assimilation on hot days and not on mild ones" begin
    hot_mean = stress_state(34.0, 0.0; steps = 0)
    hot_course = stress_state(34.0, 12.0; steps = 24)
    @test hot_course < hot_mean
    mild_mean = stress_state(20.0, 0.0; steps = 0)
    mild_course = stress_state(20.0, 12.0; steps = 24)
    # Near the optimum the curvature is small and the two nearly agree; the
    # point of the mechanism is that they diverge in the tail, not everywhere.
    @test abs(mild_course - mild_mean) < abs(hot_course - hot_mean)
    # A wider day costs more, monotonically.
    @test stress_state(34.0, 16.0; steps = 24) < stress_state(34.0, 12.0; steps = 24)
    @test stress_state(34.0, 12.0; steps = 24) < stress_state(34.0, 6.0; steps = 24)
end

@testset "the flag reaches the model and the range is demanded" begin
    @test :diurnal_temperature_stress in
          Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation)))
    @test :diurnal_temperature_stress in Base.kwarg_decl(first(methods(Agrocosm._daily_crop!)))
    @test :diurnal_temperature_stress in fieldnames(Agrocosm.SimulationConfiguration)
    configuration = Agrocosm.ablation_diurnal_stress_configuration()
    @test configuration.diurnal_temperature_stress === true
    @test configuration.subdaily_photosynthesis === false
    accepted = Set(Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation))))
    for key in keys(configuration)
        @test key in accepted
    end
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "simulations",
                           "daily_crop.jl"), String)
    @test occursin("diurnal temperature stress requires a `diurnal_range`", source)
end
