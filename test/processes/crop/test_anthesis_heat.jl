using Test
using Agrocosm

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Anthesis heat sterility on an absolute daily-maximum air temperature. The
# thresholds were measured before the mechanism was written, so the tests that
# matter are about the SHAPE - window edges, the daily-maximum reconstruction,
# inertness at rate zero, and order-independence against the three mechanisms
# that share `grain_set_fraction`.

const T32 = Float32

"""A CFT with only the named fields overridden, following the idiom the other
mechanism tests use so a rename breaks the test rather than silently passing."""
function heat_cft(; rate, threshold = 38.0, start = 0.45, stop = 0.70)
    Agrocosm.CFTParameters{T32, Int32}(;
        (f => (f === :heat_day_rate ? T32(rate) :
               f === :heat_day_temperature ? T32(threshold) :
               f === :flowering_start ? T32(start) :
               f === :flowering_end ? T32(stop) :
               getfield(Agrocosm.cft3, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
end

@testset "the loss is linear in the excess and zero below the threshold" begin
    loss(tmax, thr, inside, rate) =
        Agrocosm.anthesis_heat_loss(T32(tmax), T32(thr), inside, T32(rate))
    @test loss(40.0, 38.0, true, 0.01) ≈ T32(0.02)
    @test loss(42.0, 38.0, true, 0.01) ≈ T32(0.04)     # linear, not a count
    @test loss(38.0, 38.0, true, 0.01) == zero(T32)    # at the threshold, nothing
    @test loss(30.0, 38.0, true, 0.01) == zero(T32)    # below it, nothing
    @test loss(45.0, 38.0, false, 0.01) == zero(T32)   # outside the window, nothing
    # A day at 42 C is not a day at 38.1 C. This is the whole reason the term is
    # linear in the excess rather than counting exceedance days.
    @test loss(42.0, 38.0, true, 0.01) > loss(38.1, 38.0, true, 0.01)
end

@testset "S1/S2: rate zero is inert" begin
    # The ablation contract: `heat_day_rate = 0` must leave `grain_set_fraction`
    # exactly where it started, so the arm is bitwise the configuration below it.
    crop = test_model_state(init_crop(4, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.grain_set_fraction .= one(T32)
    state.phenology.is_growing .= Int32(1)
    Agrocosm.crop_phenology_auxiliary(crop).fphu .= T32(0.575)   # mid-window
    temperature = fill(T32(35.0), 4)
    range = fill(T32(14.0), 4)                                   # daily max 42 C
    Agrocosm.anthesis_heat!(heat_cft(rate = 0.0), crop, temperature, range)
    @test all(state.phenology.grain_set_fraction .== one(T32))

    # And with a rate it must actually move, or the zero test is vacuous.
    Agrocosm.anthesis_heat!(heat_cft(rate = 0.01), crop, temperature, range)
    @test all(state.phenology.grain_set_fraction .< one(T32))
end

@testset "the daily maximum is the mean plus half the range" begin
    # The forcing carries a mean and a range, never a maximum, so the
    # reconstruction is where a factor-of-two error would hide.
    function after(mean, range; rate = 0.01, threshold = 38.0)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(T32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= T32(0.575)
        Agrocosm.anthesis_heat!(heat_cft(rate = rate, threshold = threshold),
                                crop, fill(T32(mean), 1), fill(T32(range), 1))
        state.phenology.grain_set_fraction[1]
    end
    # mean 35 + range 14 / 2 = 42, so 4 degrees of excess at rate 0.01 -> 0.04
    @test after(35.0, 14.0) ≈ one(T32) - T32(0.04)
    # The same maximum reached two ways must give the same loss.
    @test after(35.0, 14.0) ≈ after(38.0, 8.0)
    # A negative range is data error, not a negative maximum: it is clamped.
    @test after(45.0, -10.0) ≈ after(45.0, 0.0)
end

@testset "the window is rectangular and open at both ends" begin
    function after_fphu(fphu; start = 0.45, stop = 0.70)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(T32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= T32(fphu)
        Agrocosm.anthesis_heat!(heat_cft(rate = 0.01, start = start, stop = stop),
                                crop, fill(T32(45.0), 1), fill(T32(0.0), 1))
        state.phenology.grain_set_fraction[1]
    end
    inside = after_fphu(0.575)
    @test inside < one(T32)
    # Rectangular, NOT the raised cosine: a day just inside the window loses the
    # same as a day at its centre.
    @test after_fphu(0.46) ≈ inside
    @test after_fphu(0.69) ≈ inside
    # Open at both ends, matching `flowering_weight`.
    @test after_fphu(0.45) == one(T32)
    @test after_fphu(0.70) == one(T32)
    @test after_fphu(0.20) == one(T32)
    @test after_fphu(0.90) == one(T32)
    # A degenerate window admits nothing rather than everything.
    @test after_fphu(0.5, start = 0.6, stop = 0.6) == one(T32)
    @test after_fphu(0.5, start = 0.7, stop = 0.6) == one(T32)
end

@testset "a stand that is not growing is untouched" begin
    crop = test_model_state(init_crop(2, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.grain_set_fraction .= one(T32)
    state.phenology.is_growing .= Int32(0)
    Agrocosm.crop_phenology_auxiliary(crop).fphu .= T32(0.575)
    Agrocosm.anthesis_heat!(heat_cft(rate = 0.01), crop,
                            fill(T32(45.0), 2), fill(T32(0.0), 2))
    @test all(state.phenology.grain_set_fraction .== one(T32))
end

@testset "order-independent against the heat sink" begin
    # Both subtract from `grain_set_fraction` and the state is clamped at zero,
    # so clamp(clamp(x - a) - b) == clamp(clamp(x - b) - a). Asserted rather than
    # trusted, exactly as the existing pairs are.
    function run(order)
        crop = test_model_state(init_crop(1, identity))
        state = Agrocosm.crop_prognostic(crop)
        state.phenology.grain_set_fraction .= one(T32)
        state.phenology.is_growing .= Int32(1)
        Agrocosm.crop_phenology_auxiliary(crop).fphu .= T32(0.575)
        Agrocosm.crop_stress_auxiliary(crop).heat_exposure_hours .= T32(6.0)
        cft = heat_cft(rate = 0.01)
        for step in order
            step === :anthesis && Agrocosm.anthesis_heat!(
                cft, crop, fill(T32(35.0), 1), fill(T32(14.0), 1))
            step === :sink && Agrocosm.reproductive_sink!(cft, crop)
        end
        state.phenology.grain_set_fraction[1]
    end
    @test run((:anthesis, :sink)) == run((:sink, :anthesis))
end

@testset "the flag reaches the model through the public entry" begin
    # The unit tests above call `anthesis_heat!` directly, which is exactly how a
    # missing keyword on `initialize_simulation` survived them and failed all 24
    # global jobs instead. Every process flag has to be threaded through the API
    # layer as well as the daily driver, so assert the whole chain: the keyword
    # is accepted, it lands in the config, and the ablation configuration that
    # names it produces it.
    @test :anthesis_heat in Base.kwarg_decl(
        first(methods(Agrocosm.initialize_simulation)))
    configuration = Agrocosm.ablation_anthesis_heat_configuration()
    @test configuration.anthesis_heat === true
    @test configuration.reproductive_sink === false
    # Every key this configuration sets must be a keyword the entry accepts, or
    # it fails only once a global run reaches a compute node.
    accepted = Set(Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation))))
    for key in keys(configuration)
        @test key in accepted
    end
end
