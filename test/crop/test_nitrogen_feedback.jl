using Agrocosm
using Test

@testset "Leaf-nitrogen Rubisco limitation feedback" begin
    n = CropNitrogen(Float64)   # ncleaf_min=1/58.8, ncleaf_ref=1/25, target_nc=1/30
    nlim(leaf_n, leaf_c) = Agrocosm.leaf_nitrogen_limitation(n, leaf_n, leaf_c)

    @testset "limitation from leaf N:C" begin
        # At the structural minimum N:C → fully limited (0); at/above the reference → unlimited (1).
        @test nlim((1 / 58.8) * 100.0, 100.0) ≈ 0.0
        @test nlim((1 / 25) * 100.0, 100.0) ≈ 1.0
        @test nlim((1 / 14.3) * 100.0, 100.0) == 1.0     # above reference → clamped to 1
        # midway between min and ref
        nc_mid = 0.5 * (1 / 58.8 + 1 / 25)
        @test nlim(nc_mid * 100.0, 100.0) ≈ 0.5 rtol = 1e-6
        # monotone increasing in leaf nitrogen
        @test nlim(2.0, 100.0) > nlim(1.0, 100.0)
    end

    @testset "bootstrap: no leaf carbon → unlimited" begin
        @test nlim(0.0, 0.0) == 1.0                      # avoids the early-growth deadlock
    end

    @testset "photosynthesis honours the limitation factor" begin
        # With nitrogen_limitation < 1, net assimilation is reduced vs the unlimited case.
        photo = CropPhotosynthesis(Float64)
        args = (12.0, 22.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)   # cmass, T, sw, pres, co2, LAI, λc, β
        _, An_full = Agrocosm.compute_respiration_assimilation(photo, args...)          # default nlim = 1
        _, An_half = Agrocosm.compute_respiration_assimilation(photo, args..., 0.5)     # nlim = 0.5
        @test An_half < An_full
        @test An_full > 0
    end
end
