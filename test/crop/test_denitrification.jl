using Agrocosm
using Test

@testset "Crop soil denitrification" begin
    dn = CropDenitrification(Float64)   # CDN=1.2, n2o_fraction=0.11

    @testset "temperature factor" begin
        @test Agrocosm.denitrification_temperature_factor(dn, 50.0) == 0.0     # above 45.9 °C
        @test Agrocosm.denitrification_temperature_factor(dn, -5.0) ≈ 0.0326   # at/below 0 °C
        # closed form for 0 < T ≤ 45.9
        t = 20.0
        expected = 0.0326 + 0.00351 * t^1.652 - (t / 41.748)^7.19
        @test Agrocosm.denitrification_temperature_factor(dn, t) ≈ expected rtol = 1e-9
    end

    @testset "moisture factor" begin
        @test Agrocosm.denitrification_moisture_factor(dn, 1.0) ≤ 1.0          # capped
        @test Agrocosm.denitrification_moisture_factor(dn, 0.8) ≈ min(1.0, 6.664096e-10 * exp(20.92912 * 0.8)) rtol = 1e-9
        # increases with wetness (until the cap)
        @test Agrocosm.denitrification_moisture_factor(dn, 0.9) > Agrocosm.denitrification_moisture_factor(dn, 0.5)
    end

    @testset "gross denitrification and N2O/N2 split" begin
        g, n2o, n2 = Agrocosm.gross_denitrification(dn, 100.0, 20.0, 0.8, 1.0)
        @test 0.0 ≤ g ≤ 100.0
        @test n2o + n2 ≈ g                    # gaseous nitrogen conserved
        @test n2o ≈ 0.11 * g
        # capped at the nitrate stock even with maximal factors
        gcap, _, _ = Agrocosm.gross_denitrification(dn, 1.0, 20.0, 1.0, 100.0)
        @test gcap ≤ 1.0
        # more organic carbon → more denitrification (up to the cap)
        g_low, _, _ = Agrocosm.gross_denitrification(dn, 100.0, 20.0, 0.8, 0.5)
        @test g > g_low
    end
end
