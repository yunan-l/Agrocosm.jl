using Test
using Agrocosm

# FAO-56's readily-available-water plateau on the transpiration supply.
#
# The contract that matters is BITWISE equivalence at `depletion_fraction = 0`:
# this changes the water supply of every simulated cell, so an accidental drift
# at zero would silently move every result in the project, including the
# baselines `docs/21` and `docs/22` are written against.

const TW = Float32

@testset "at p = 0 the supply is bitwise the LPJmL form" begin
    # The expression this replaces, written out, so a refactor that changes the
    # arithmetic fails here rather than moving every yield in the project.
    lpjml(emax, wr, rc) = emax * wr * (one(TW) - exp(TW(-0.04) * rc))
    for wr in TW[0.0, 0.013, 0.25, 0.5, 0.732, 0.999, 1.0],
        rc in TW[0.0, 3.7, 50.0, 300.0],
        emax in TW[1.0, 5.0, 7.0]
        @test Agrocosm.compute_transpiration_supply(emax, wr, rc, zero(TW)) ===
              lpjml(emax, wr, rc)
        # The default argument must also be zero, or existing callers drift.
        @test Agrocosm.compute_transpiration_supply(emax, wr, rc) ===
              lpjml(emax, wr, rc)
    end
    # And the fraction itself is the identity at zero, for every input.
    for wr in TW[0.0, 0.1, 0.37, 0.5, 0.9, 1.0]
        @test Agrocosm.compute_available_fraction(wr, zero(TW)) === wr
    end
end

@testset "the plateau is FAO-56's readily available water" begin
    f(wr, p) = Agrocosm.compute_available_fraction(TW(wr), TW(p))
    # With p = 0.5 the crop is unstressed until half the available water is gone.
    @test f(1.0, 0.5) == one(TW)
    @test f(0.6, 0.5) == one(TW)
    @test f(0.5, 0.5) == one(TW)     # exactly at the threshold, still unstressed
    @test f(0.25, 0.5) ≈ TW(0.5)     # halfway through the remainder
    @test f(0.0, 0.5) == zero(TW)
    # Larger p means a longer plateau, so never less supply at the same water.
    for wr in TW[0.05, 0.2, 0.4, 0.7]
        @test f(wr, 0.55) >= f(wr, 0.2) >= f(wr, 0.0)
    end
    # Monotone in water, and clamped at one.
    @test f(0.3, 0.55) < f(0.4, 0.55) < f(0.5, 0.55)
    @test f(1.0, 0.55) == one(TW)
    # A degenerate p = 1 means "never stressed while there is any water", not Inf.
    @test f(0.1, 1.0) == one(TW)
    @test f(0.0, 1.0) == zero(TW)
    @test isfinite(f(0.1, 1.0))
end

@testset "the published values are FAO-56 Table 22" begin
    @test Agrocosm.fao56_depletion_fraction(1) == 0.55   # wheat
    @test Agrocosm.fao56_depletion_fraction(2) == 0.20   # rice, lowest in the table
    @test Agrocosm.fao56_depletion_fraction(3) == 0.55   # maize, field grain
    @test Agrocosm.fao56_depletion_fraction(9) == 0.50   # soybean
    # A CFT the table does not cover leaves the supply unmodified rather than
    # guessing a value.
    @test Agrocosm.fao56_depletion_fraction(12) == 0.0
    # Every shipped CFT is INERT: turning the plateau on has to be an explicit
    # experiment, because it changes every simulated yield.
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.depletion_fraction == zero(TW)
    end
end

@testset "the parameter reaches the model and relieves stress" begin
    # A function test is not enough: `depletion_fraction` has to travel from the
    # CFT into both kernels that compute supply. Asserting on the field alone
    # would pass with the parameter unused.
    @test :depletion_fraction in fieldnames(Agrocosm.CFTParameters)
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "processes", "crop",
                           "transpiration.jl"), String)
    # Both call sites, and the LAI senescence scalar, must carry it - `wscal`
    # drives senescence and must not use a different water stress from
    # allocation.
    @test count("depletion_fraction", source) >= 6
    @test occursin("compute_available_fraction(wr, T(depletion_fraction))", source)
end
