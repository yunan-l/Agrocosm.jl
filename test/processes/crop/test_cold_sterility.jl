using Test
using Agrocosm

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Cold sterility on an absolute daily-minimum air temperature: the mirror of
# `anthesis_heat!`. The tests that matter are about SHAPE - the sign of the
# deficit, the daily-minimum reconstruction, its own window rather than the
# flowering one, inertness at rate zero, and order-independence against the
# mechanisms that share `grain_set_fraction`.

const TC32 = Float32

"""A CFT with only the named fields overridden, so a rename breaks the test."""
function cold_cft(; rate, threshold = 17.0, start = 0.30, stop = 0.70)
    Agrocosm.CFTParameters{TC32, Int32}(;
        (f => (f === :cold_night_rate ? TC32(rate) :
               f === :cold_night_temperature ? TC32(threshold) :
               f === :cold_start ? TC32(start) :
               f === :cold_end ? TC32(stop) :
               getfield(Agrocosm.cft2, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
end

@testset "the loss is linear in the deficit and zero above the threshold" begin
    loss(tmin, thr, inside, rate) =
        Agrocosm.cold_sterility_loss(TC32(tmin), TC32(thr), inside, TC32(rate))
    @test loss(15.0, 17.0, true, 0.01) ≈ TC32(0.02)
    @test loss(13.0, 17.0, true, 0.01) ≈ TC32(0.04)    # linear, not a count
    @test loss(17.0, 17.0, true, 0.01) == zero(TC32)   # at the threshold, nothing
    @test loss(25.0, 17.0, true, 0.01) == zero(TC32)   # above it, nothing
    @test loss(5.0, 17.0, false, 0.01) == zero(TC32)   # outside the window, nothing
    # The sign is the whole difference from the heat term, so assert it rather
    # than trust it: COLDER must cost more, not less.
    @test loss(10.0, 17.0, true, 0.01) > loss(16.9, 17.0, true, 0.01)
    # And a sub-zero minimum is an ordinary deficit, not a special case.
    @test loss(-3.0, 0.0, true, 0.01) ≈ TC32(0.03)
end

@testset "S1/S2: rate zero is inert" begin
    # The ablation contract: `cold_night_rate = 0` must leave the state exactly
    # where it started, so the arm is bitwise the configuration below it.
    crop = test_model_state(init_crop(4, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.grain_set_fraction .= one(TC32)
    state.phenology.is_growing .= Int32(1)
    Agrocosm.crop_phenology_auxiliary(crop).fphu .= TC32(0.5)
    temperature = fill(TC32(19.0), 4)
    range = fill(TC32(12.0), 4)                        # daily min 13 C
    Agrocosm.cold_sterility!(cold_cft(rate = 0.0), crop, temperature, range)
    @test all(state.phenology.grain_set_fraction .== one(TC32))

    # And with a rate it must actually move, or the zero test is vacuous.
    Agrocosm.cold_sterility!(cold_cft(rate = 0.01), crop, temperature, range)
    @test all(state.phenology.grain_set_fraction .< one(TC32))
end

@testset "the daily minimum is the mean MINUS half the range" begin
    # The sign of this reconstruction is where the copy-paste from the heat
    # kernel would break, and it would break silently: a mean of 19 with a range
    # of 12 is a minimum of 13, not 25.
    function after(mean, range; rate = 0.01, threshold = 17.0)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(TC32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= TC32(0.5)
        Agrocosm.cold_sterility!(cold_cft(rate = rate, threshold = threshold),
                                 crop, fill(TC32(mean), 1), fill(TC32(range), 1))
        state.phenology.grain_set_fraction[1]
    end
    # mean 19 - range 12 / 2 = 13, so 4 degrees of deficit at 0.01 -> 0.04
    @test after(19.0, 12.0) ≈ one(TC32) - TC32(0.04)
    # The same minimum reached two ways must give the same loss.
    @test after(19.0, 12.0) ≈ after(16.0, 6.0)
    # WIDENING the range must cost MORE, which is the opposite of the heat term
    # and the assertion a sign error would fail.
    @test after(19.0, 16.0) < after(19.0, 12.0)
    # A negative range is data error, not a positive minimum: it is clamped.
    @test after(10.0, -10.0) ≈ after(10.0, 0.0)
end

@testset "the cold window is its own, not the flowering window" begin
    function after_fphu(fphu; start = 0.30, stop = 0.70)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(TC32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= TC32(fphu)
        Agrocosm.cold_sterility!(cold_cft(rate = 0.01, start = start, stop = stop),
                                 crop, fill(TC32(5.0), 1), fill(TC32(0.0), 1))
        state.phenology.grain_set_fraction[1]
    end
    inside = after_fphu(0.5)
    @test inside < one(TC32)
    # Rectangular: a day just inside costs the same as one at the centre.
    @test after_fphu(0.31) ≈ inside
    @test after_fphu(0.69) ≈ inside
    # Open at both ends, matching `anthesis_heat!`.
    @test after_fphu(0.30) == one(TC32)
    @test after_fphu(0.70) == one(TC32)
    @test after_fphu(0.90) == one(TC32)
    # The point of a separate window: microsporogenesis precedes anthesis, so
    # fphu 0.35 is cold-sensitive and NOT heat-sensitive. If the kernel ever
    # reads `flowering_start` instead, this is the assertion that fails.
    @test after_fphu(0.35) < one(TC32)
    @test Agrocosm.cft2.cold_start < Agrocosm.cft2.flowering_start
    # A degenerate window admits nothing rather than everything.
    @test after_fphu(0.5, start = 0.6, stop = 0.6) == one(TC32)
    @test after_fphu(0.5, start = 0.7, stop = 0.6) == one(TC32)
end

@testset "a stand that is not growing is untouched" begin
    crop = test_model_state(init_crop(2, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.grain_set_fraction .= one(TC32)
    state.phenology.is_growing .= Int32(0)
    Agrocosm.crop_phenology_auxiliary(crop).fphu .= TC32(0.5)
    Agrocosm.cold_sterility!(cold_cft(rate = 0.01), crop,
                             fill(TC32(0.0), 2), fill(TC32(0.0), 2))
    @test all(state.phenology.grain_set_fraction .== one(TC32))
end

@testset "order-independent against the heat sink and anthesis heat" begin
    # All three subtract from `grain_set_fraction` and the state is clamped, so
    # the composition is order-free. Asserted rather than trusted.
    function run(order)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(TC32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= TC32(0.5)
        Agrocosm.crop_stress_auxiliary(crop).heat_exposure_hours .= TC32(6.0)
        cft = cold_cft(rate = 0.01)
        for step in order
            step === :cold && Agrocosm.cold_sterility!(
                cft, crop, fill(TC32(19.0), 1), fill(TC32(12.0), 1))
            step === :sink && Agrocosm.reproductive_sink!(cft, crop)
        end
        state.phenology.grain_set_fraction[1]
    end
    # To an ulp, not bitwise: float subtraction is not associative, so
    # `(x - a) - b` and `(x - b) - a` may differ in the last bit. The contract is
    # that the CLAMP composes - no order produces a different amount of damage -
    # not that the arithmetic is bit-identical.
    @test run((:cold, :sink)) ≈ run((:sink, :cold)) atol = eps(TC32)
end

@testset "the per-crop thresholds are the ones the sweep and the literature give" begin
    # These are the numbers `docs/19` records, and a silent edit to `cft.jl`
    # should break a test rather than change a published mechanism.
    @test Agrocosm.cft1.cold_night_temperature == TC32(0.0)    # wheat, frost at anthesis
    @test Agrocosm.cft2.cold_night_temperature == TC32(17.0)   # rice, LITERATURE not sweep
    @test Agrocosm.cft3.cold_night_temperature == TC32(2.0)    # maize, strongest cold signal
    @test Agrocosm.cft9.cold_night_temperature == TC32(0.0)    # soybean, weakest
    # Ships AT its bound so `rate_scale` can scale it down; the mechanism is kept
    # out of production runs by the process flag, not by a zero rate.
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.cold_night_rate == TC32(0.04)
        @test cft.cold_start == TC32(0.30)
        @test cft.cold_end == TC32(0.70)
    end
end

@testset "the flag reaches the model through the public entry" begin
    # A missing keyword on `initialize_simulation` once survived every unit test
    # here and failed all 24 global jobs instead, because the units call the
    # kernel directly. Assert the whole chain.
    @test :cold_sterility in Base.kwarg_decl(
        first(methods(Agrocosm.initialize_simulation)))
    configuration = Agrocosm.ablation_cold_sterility_configuration()
    @test configuration.cold_sterility === true
    @test configuration.reproductive_sink === false
    accepted = Set(Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation))))
    for key in keys(configuration)
        @test key in accepted
    end
    # It is a fixed comparison cell, so it must refuse to have its own flag set.
    @test_throws ArgumentError Agrocosm.ablation_cold_sterility_configuration(
        cold_sterility = false)
    @test_throws ArgumentError Agrocosm.ablation_cold_sterility_configuration(
        organ_temperature = true)
end

@testset "the daily driver accepts the flag and the runtime records it" begin
    # `daily_crop_C3!` forwards `kwargs...`, so a keyword test there is vacuous -
    # the driver that actually declares the flag is `_daily_crop!`, and the
    # metadata struct is what a run is reproduced from.
    @test :cold_sterility in Base.kwarg_decl(first(methods(Agrocosm._daily_crop!)))
    # The configuration is what a checkpoint is validated against and what the
    # run manifest records, so the flag has to survive into it.
    @test :cold_sterility in fieldnames(Agrocosm.SimulationConfiguration)
end
