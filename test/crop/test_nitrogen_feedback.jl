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

    @testset "photosynthesis honours the nitrogen capacity" begin
        # Capping the Rubisco capacity below the light-derived Vc_max reduces net assimilation.
        photo = CropPhotosynthesis(Float64)
        args = (12.0, 22.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)   # cmass, T, sw, pres, co2, LAI, λc, β
        _, An_full = Agrocosm.compute_respiration_assimilation(photo, args...)          # default capacity = Inf (no cap)
        _, An_capped = Agrocosm.compute_respiration_assimilation(photo, args..., 0.0)   # zero capacity → Vc_max capped to 0
        @test An_capped < An_full
        @test An_full > 0
    end

    @testset "nitrogen-supported Vc_max capacity" begin
        limit = CropNitrogenVcmaxLimit(Float64)
        # More leaf nitrogen (at fixed leaf carbon) supports a larger Rubisco capacity; with no leaf
        # carbon yet the capacity is Inf, so early growth is not nitrogen-deadlocked.
        cap(leaf_n, leaf_c) = Agrocosm.nitrogen_supported_vcmax(limit, leaf_n, leaf_c, 25.0)
        @test cap(0.003, 0.1) > cap(0.002, 0.1) > 0
        @test cap(0.003, 0.0) == Inf
    end
end
