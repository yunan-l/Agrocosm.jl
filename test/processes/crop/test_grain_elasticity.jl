using Test
using Agrocosm
using Agrocosm: published_grain_traits, published_grain_calibration,
                grain_traits_from_elasticity, saturating_grain_number

@testset "the shipped traits imply exactly one elasticity" begin
    # The shipped pair was derived from a LEVEL constraint alone. Re-deriving it
    # from the level plus the elasticity it implicitly carries must return it
    # unchanged - that is what makes the elasticity the free parameter it was.
    for cft_id in (1, 2, 3, 9)
        half, per_grain, ceiling = published_grain_traits(cft_id)
        grains, npp = published_grain_calibration(cft_id)
        implied = half / (npp + half)
        derived = grain_traits_from_elasticity(cft_id, implied)
        # To three significant figures: the shipped `half` values carry one
        # decimal, so the round trip cannot be tighter than their own rounding
        # (rice returns 39988.6 against a tabled 40000).
        @test derived[1] ≈ half rtol = 1e-3
        @test derived[2] === per_grain
        @test derived[3] ≈ ceiling rtol = 1e-3
        # And the level constraint itself: typical window NPP returns the
        # typical published grain number, whatever the elasticity.
        for e in (implied, 0.4, 0.53, 0.7)
            h, _, c = grain_traits_from_elasticity(cft_id, e)
            @test saturating_grain_number(npp, c, h) ≈ grains rtol = 1e-6
        end
        # The implied elasticities are NOT a shared constant: 0.25 for wheat and
        # rice, 0.29 for soybean, 0.33 for maize. Nothing chose them.
        @test 0.24 < implied < 0.34
    end
end

@testset "a higher elasticity costs more grain for the same shortfall" begin
    # Braunschweig maize 2008: window assimilate fell to 0.543 of the wet arm and
    # grain number to 0.723. The shipped curve turns that shortfall into 0.806.
    _, npp = published_grain_calibration(3)
    shortfall = 0.543
    function ratio(e)
        h, _, c = grain_traits_from_elasticity(3, e)
        saturating_grain_number(npp * shortfall, c, h) / saturating_grain_number(npp, c, h)
    end
    # Evaluated at the TYPICAL window NPP, so these differ a little from the run
    # itself (0.806 and 0.714), where the dry arm's operating point sits lower on
    # the curve and its local elasticity is higher.
    @test ratio(1 / 3) > 0.77          # the shipped curve is too forgiving
    @test 0.67 < ratio(0.53) < 0.72    # the measured one lands near the measurement
    @test ratio(0.7) < ratio(0.53) < ratio(1 / 3)
end

@testset "the derivation refuses what it cannot do" begin
    # An elasticity outside (0,1) has no saturating curve, and a crop the table
    # does not cover returns a zero ceiling, which switches the sink off rather
    # than guessing.
    @test grain_traits_from_elasticity(3, 0.0)[3] == 0.0
    @test grain_traits_from_elasticity(3, 1.0)[3] == 0.0
    @test grain_traits_from_elasticity(7, 0.53) == (0.0, 0.0, 0.0)
    @test published_grain_calibration(7) == (0.0, 0.0)
end
