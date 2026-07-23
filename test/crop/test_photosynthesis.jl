using Agrocosm
using Terrarium
using Test

# Shorthand for the scalar primitive.
ra(photo, T, sw, pres, co2, LAI, λc, β; cmass = 12.0) =
    Agrocosm.compute_respiration_assimilation(photo, cmass, T, sw, pres, co2, LAI, λc, β)

@testset "Crop photosynthesis" begin
    @testset "C3 reproduces Terrarium LUEPhotosynthesis" begin
        # With PALADYN/BIOME3 defaults, CropPhotosynthesis C3 is the same biochemistry as
        # Terrarium's LUEPhotosynthesis. Cross-check the scalar (Rd, An) across an input grid.
        crop = CropPhotosynthesis(Float64)               # C3, LUE defaults
        lue = LUEPhotosynthesis(Float64)
        mat = MaterialConstants(Float64)
        for T in (0.0, 10.0, 20.0, 30.0), sw in (50.0, 200.0, 600.0),
                LAI in (0.5, 2.0, 5.0), λc in (0.6, 0.7, 0.8), β in (0.3, 1.0)
            pres, co2 = 1.0e5, 380.0
            Rd_c, An_c = ra(crop, T, sw, pres, co2, LAI, λc, β)
            Rd_l, An_l = Terrarium.compute_respiration_assimilation(lue, mat, T, sw, pres, co2, LAI, λc, β)
            @test Rd_c ≈ Rd_l rtol = 1e-10
            @test An_c ≈ An_l rtol = 1e-10
        end
    end

    @testset "gating invariants (C3 and C4)" begin
        for path in (C3Pathway(), C4Pathway())
            photo = CropPhotosynthesis(Float64; pathway = path)
            # No light → no assimilation or respiration.
            @test ra(photo, 20.0, 0.0, 1.0e5, 380.0, 3.0, 0.7, 1.0) == (0.0, 0.0)
            # No leaves → nothing.
            @test ra(photo, 20.0, 400.0, 1.0e5, 380.0, 0.0, 0.7, 1.0) == (0.0, 0.0)
            # Too cold → nothing.
            @test ra(photo, -5.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0) == (0.0, 0.0)
            # Active conditions → net assimilation is finite and non-negative.
            Rd, An = ra(photo, 22.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)
            @test isfinite(Rd) && isfinite(An)
            @test An ≥ 0
        end
    end

    @testset "C4 behaviour" begin
        c4 = CropPhotosynthesis(Float64; pathway = C4Pathway())
        # Net assimilation increases with incoming light.
        _, An_low = ra(c4, 25.0, 100.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)
        _, An_high = ra(c4, 25.0, 600.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)
        @test An_high > An_low
        # φ = min(1, λc/λ_mc4) saturates at λc ≥ λ_mc4 = 0.4: gross terms match for λc = 0.4 and 0.8.
        # (compare c₁ via the light-limited branch at low LAI where JE dominates)
        Rd1, An1 = ra(c4, 25.0, 300.0, 1.0e5, 380.0, 3.0, 0.5, 1.0)
        Rd2, An2 = ra(c4, 25.0, 300.0, 1.0e5, 380.0, 3.0, 0.9, 1.0)
        @test An1 ≈ An2 rtol = 1e-12   # saturated for both λc ≥ 0.4
        # C4 high-temperature cutoff is 55 °C (vs C3's 45 °C): active at 50 °C.
        _, An_hot = ra(c4, 50.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0)
        @test An_hot ≥ 0
        c3 = CropPhotosynthesis(Float64; pathway = C3Pathway())
        @test ra(c3, 50.0, 400.0, 1.0e5, 380.0, 3.0, 0.7, 1.0) == (0.0, 0.0)   # C3 cut off at 45
    end
end
