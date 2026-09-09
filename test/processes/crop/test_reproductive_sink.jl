using Agrocosm
using Test

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Gates S1-S7 of docs/05_reproductive_sink_design.md. S8 (Enzyme) lives with the
# AD suite in test/ad/test_weather_attribution.jl.

@testset "Flowering window weight" begin
    start, stop = 0.45, 0.70
    # Outside the window heat does nothing to grain set: the sensitivity is
    # developmental, which is the whole point of having a window.
    for fphu in (0.0, 0.2, 0.44, 0.45, 0.70, 0.85, 1.0)
        @test Agrocosm.flowering_weight(fphu, start, stop) == 0.0
    end
    for fphu in (0.46, 0.5, 0.575, 0.65, 0.69)
        @test Agrocosm.flowering_weight(fphu, start, stop) > 0
    end
    # Peaks mid-window and is symmetric about it.
    peak = Agrocosm.flowering_weight(0.575, start, stop)
    @test peak ≈ 1.0
    for offset in (0.02, 0.05, 0.1)
        @test Agrocosm.flowering_weight(0.575 - offset, start, stop) ≈
              Agrocosm.flowering_weight(0.575 + offset, start, stop)
    end
    # Continuous at both edges: a hard window would put kinks in the reverse pass.
    @test Agrocosm.flowering_weight(0.4501, start, stop) < 1e-3
    @test Agrocosm.flowering_weight(0.6999, start, stop) < 1e-3
    # Degenerate window is inert rather than an error.
    @test Agrocosm.flowering_weight(0.5, 0.7, 0.7) == 0.0
    @test Agrocosm.flowering_weight(0.5, 0.8, 0.6) == 0.0
end

@testset "Grain-set loss is linear and non-negative" begin
    @test Agrocosm.grain_set_loss(0.0, 1.0, 0.02) == 0.0
    @test Agrocosm.grain_set_loss(5.0, 0.0, 0.02) == 0.0
    @test Agrocosm.grain_set_loss(5.0, 1.0, 0.0) == 0.0
    @test Agrocosm.grain_set_loss(4.0, 1.0, 0.02) ≈ 0.08
    # Linear in exposure, so doubling the hours doubles the loss. This is what
    # makes duration, rather than peak, the currency.
    @test Agrocosm.grain_set_loss(8.0, 1.0, 0.02) ≈
          2 * Agrocosm.grain_set_loss(4.0, 1.0, 0.02)
    @test Agrocosm.grain_set_loss(-1.0, 1.0, 0.02) == 0.0
end

"""Step the sink kernel over a sequence of days and return the trajectory."""
function sterility_trajectory(exposures, fphus; rate = 0.02, growing = true,
                              start = 0.45, stop = 0.70)
    T = Float32
    cft = Agrocosm.CFTParameters{T, Int32}(;
        (f => (f === :sterility_rate ? T(rate) :
               f === :flowering_start ? T(start) :
               f === :flowering_end ? T(stop) :
               getfield(Agrocosm.cft1, f)) for f in fieldnames(Agrocosm.CFTParameters))...)
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    Agrocosm.crop_prognostic(state).phenology.grain_set_fraction .= one(T)
    Agrocosm.crop_prognostic(state).phenology.is_growing .= Int32(growing)
    trajectory = T[]
    for (exposure, fphu) in zip(exposures, fphus)
        Agrocosm.crop_stress_auxiliary(state).heat_exposure_hours .= T(exposure)
        Agrocosm.crop_phenology_auxiliary(state).fphu .= T(fphu)
        Agrocosm.reproductive_sink!(cft, state)
        push!(trajectory, Agrocosm.crop_prognostic(state).phenology.grain_set_fraction[1])
    end
    return trajectory
end

@testset "S4: sterility is monotone and irreversible" begin
    # A heat wave during flowering, then a cool spell. Grain that failed to set
    # is not recovered - unlike every other stress path in the model, which
    # relaxes as soon as the stress lifts.
    hot_then_cool = sterility_trajectory(
        [6.0, 6.0, 6.0, 0.0, 0.0, 0.0], fill(0.575, 6),
    )
    @test issorted(hot_then_cool; rev = true)
    @test hot_then_cool[3] < 1.0
    @test hot_then_cool[end] == hot_then_cool[3]   # cool days recover nothing

    # Never leaves [0, 1], even under absurd exposure.
    crushed = sterility_trajectory(fill(500.0, 5), fill(0.575, 5))
    @test all(0.0 .<= crushed .<= 1.0)
    @test crushed[end] == 0.0
    @test issorted(crushed; rev = true)
end

@testset "S5: sensitivity is developmental, not seasonal" begin
    # Identical heat, inside and outside the flowering window.
    inside = sterility_trajectory(fill(6.0, 4), fill(0.575, 4))
    before = sterility_trajectory(fill(6.0, 4), fill(0.20, 4))
    after = sterility_trajectory(fill(6.0, 4), fill(0.90, 4))
    @test inside[end] < 1.0
    @test before[end] == 1.0
    @test after[end] == 1.0

    # And nothing happens at all when no crop is standing.
    @test sterility_trajectory(fill(6.0, 4), fill(0.575, 4); growing = false)[end] == 1.0
end

@testset "S6: duration matters, not the daily peak" begin
    # Two days that a daily-maximum criterion could not tell apart: same
    # threshold exceeded, different hours above it. This is the Fig 1
    # temporal-aggregation argument applied to the sink.
    brief = sterility_trajectory([1.0], [0.575])
    sustained = sterility_trajectory([8.0], [0.575])
    @test sustained[1] < brief[1] < 1.0
    @test (1 - sustained[1]) ≈ 8 * (1 - brief[1]) rtol = 1e-5
end

@testset "S1/S2: off and zero-rate are inert" begin
    # Rate zero carries the state variable without ever moving it: the exact
    # "sink off" arm of the ablation experiment.
    @test sterility_trajectory(fill(10.0, 5), fill(0.575, 5); rate = 0.0) == ones(Float32, 5)
end

@testset "Grain set is initialised to one" begin
    # The state has to start at full set, not at the zero every other float
    # field is initialised to, or the first season would be born sterile.
    for T in (Float32, Float64)
        phenology = Agrocosm.init_crop_phenology(T, 4, identity)
        @test all(phenology.grain_set_fraction .== one(T))
    end
end

# Reset at sowing is a multi-season property and cannot be seen in a single
# kernel call, so it is verified end to end rather than here: without it the
# second season inherits the first season's sterility. See the multi-season
# assertion in the organ-temperature/sink end-to-end run.
