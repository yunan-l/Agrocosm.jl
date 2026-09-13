using Test
using Agrocosm

# The CERES-style grain sink that replaces the prescribed harvest index.
#
# `docs/34` measured why: the inherited index responds only to season water
# sufficiency, through a logistic 97.8% saturated at the lowest value real cells
# reach, so across the whole realized range it moves through 2% of its own span
# while BINDING 83% of maize days. Yield was a fixed fraction of biomass, which
# makes every yield loss that is not a biomass loss structurally unrepresentable.

const TG = Float32

@testset "the mechanism ships inert and is bitwise the inherited index" begin
    # `grain_number_half_carbon = 0` takes the old branch, unchanged.
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.grain_number_half_carbon == zero(TG)
    end
    @test Agrocosm.saturating_grain_number(TG(500), TG(4500), zero(TG)) == zero(TG)
end

@testset "grain number saturates toward the ceiling and never exceeds it" begin
    ceiling, half = TG(4500), TG(86.1)
    n(x) = Agrocosm.saturating_grain_number(TG(x), ceiling, half)
    @test n(0) == zero(TG)
    @test n(half) ≈ ceiling / 2
    # Monotone, bounded, and still moving at the top - the property a hard cap
    # lacked. Measured before this form was adopted: linear-plus-cap pinned rice,
    # maize and irrigated wheat at exactly `ceiling * maximum_grain_carbon` every
    # season, which is a constant yield and the defect this replaces.
    @test n(10) < n(100) < n(1000) < n(10_000) < ceiling
    # Bounded BY the ceiling, not strictly below it: at 1e9 gC the ratio rounds
    # to exactly 1 in Float32, so `n` returns the ceiling itself. That is the
    # correct behaviour for an asymptote in finite precision, and asserting a
    # strict inequality there was asserting something untrue.
    @test n(1e9) <= ceiling
    @test n(1e30) <= ceiling
    # Negative assimilate is not a negative grain number.
    @test n(-50) == zero(TG)
end

@testset "the published traits reproduce the published grain numbers" begin
    # Each crop's `half_carbon` is fixed by requiring the TYPICAL window NPP this
    # model produces to return the TYPICAL published grain number. No yield was
    # consulted, so this asserts a derivation rather than a fit.
    for (cft_id, typical_npp, typical_grains) in (
            (1, 126.9, 15000), (2, 88.0, 30000), (3, 172.2, 3000), (9, 153.4, 2500))
        half, per_grain, ceiling = Agrocosm.published_grain_traits(cft_id)
        got = Agrocosm.saturating_grain_number(TG(typical_npp), TG(ceiling), TG(half))
        @test isapprox(got, TG(typical_grains); rtol = 1e-3)
        @test ceiling > typical_grains          # a ceiling, not an operating point
        @test 0 < per_grain < 1                 # gC per grain, not milligrams
    end
    @test Agrocosm.published_grain_traits(99) == (0.0, 0.0, 0.0)
end

@testset "the sink is grains times what each can still hold" begin
    sink(n, w, fphu, stop, fill) =
        Agrocosm.grain_sink_carbon(TG(n), TG(w), TG(fphu), TG(stop), TG(fill))
    # Nothing before the window closes, everything once filling completes.
    @test sink(3000, 0.126, 0.70, 0.70, 1.0) == zero(TG)
    @test sink(3000, 0.126, 1.00, 0.70, 1.0) ≈ TG(3000 * 0.126)
    @test sink(3000, 0.126, 0.85, 0.70, 1.0) ≈ TG(3000 * 0.126) / 2
    # Linear in grain number, so a lost grain is a lost grain.
    @test sink(1500, 0.126, 1.0, 0.70, 1.0) ≈ sink(3000, 0.126, 1.0, 0.70, 1.0) / 2
    # `grain_fill_fraction` scales the WEIGHT here; `grain_set_fraction` scales
    # the NUMBER at the call site. That separation is the whole point: until now
    # both multiplied the same constant index, so neither acted on what it names.
    @test sink(3000, 0.126, 1.0, 0.70, 0.5) ≈ sink(3000, 0.126, 1.0, 0.70, 1.0) / 2
    # A degenerate window does not divide by zero.
    @test sink(3000, 0.126, 1.0, 1.0, 1.0) ≈ TG(3000 * 0.126)
end
