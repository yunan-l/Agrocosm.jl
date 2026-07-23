using Agrocosm
using Test

@testset "Crop stomatal conductance" begin
    sc = CropStomatalConductance(Float64)   # λ_opt=0.8, λ_min=0.2, g_min=1e-4, g_max=1e-2, k_ext=0.5

    @testset "leaf-to-air CO₂ ratio vs water stress" begin
        # Well-watered → LAMBDA_OPT; fully stressed → λ_min; linear and monotone in β.
        @test Agrocosm.compute_leaf_to_air_co2_ratio(sc, 1.0) ≈ 0.8
        @test Agrocosm.compute_leaf_to_air_co2_ratio(sc, 0.0) ≈ 0.2
        @test Agrocosm.compute_leaf_to_air_co2_ratio(sc, 0.5) ≈ 0.5
        λ = [Agrocosm.compute_leaf_to_air_co2_ratio(sc, β) for β in 0.0:0.1:1.0]
        @test issorted(λ)
        @test all(0.2 .≤ λ .≤ 0.8)
    end

    @testset "canopy conductance bounds and monotonicity" begin
        # No canopy or no water → minimum conductance.
        @test Agrocosm.compute_canopy_conductance(sc, 0.0, 1.0) ≈ 1.0e-4   # LAI = 0
        @test Agrocosm.compute_canopy_conductance(sc, 3.0, 0.0) ≈ 1.0e-4   # β = 0
        # Increases with LAI and with β, stays within [g_min, g_max].
        @test Agrocosm.compute_canopy_conductance(sc, 5.0, 1.0) > Agrocosm.compute_canopy_conductance(sc, 1.0, 1.0)
        @test Agrocosm.compute_canopy_conductance(sc, 3.0, 1.0) > Agrocosm.compute_canopy_conductance(sc, 3.0, 0.5)
        for LAI in (0.0, 1.0, 3.0, 8.0), β in (0.0, 0.5, 1.0)
            g = Agrocosm.compute_canopy_conductance(sc, LAI, β)
            @test 1.0e-4 ≤ g ≤ 1.0e-2
        end
    end
end
