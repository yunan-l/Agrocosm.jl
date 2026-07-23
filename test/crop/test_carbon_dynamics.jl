using Agrocosm
using Test

@testset "Crop carbon dynamics" begin
    cd = CropCarbonDynamics(Float64)   # SLA=10, awl=2, LAI_min=1, LAI_max=6, γ...

    @testset "balanced LAI from carbon pool" begin
        # LAI_b = C_veg / (2/SLA + awl) = C_veg / 2.2, clamped ≥ 0.
        @test Agrocosm.compute_balanced_lai(cd, 0.0) == 0.0
        @test Agrocosm.compute_balanced_lai(cd, 0.5) ≈ 0.5 / 2.2
        @test Agrocosm.compute_balanced_lai(cd, 2.2) ≈ 1.0
        @test Agrocosm.compute_balanced_lai(cd, -1.0) == 0.0        # defensive clamp
        # monotone increasing in carbon
        @test Agrocosm.compute_balanced_lai(cd, 3.0) > Agrocosm.compute_balanced_lai(cd, 1.0)
    end

    @testset "NPP partitioning factor" begin
        @test Agrocosm.compute_lambda_npp(cd, 0.5) == 0.0          # below LAI_min
        @test Agrocosm.compute_lambda_npp(cd, 1.0) ≈ 0.0          # at LAI_min
        @test Agrocosm.compute_lambda_npp(cd, 3.5) ≈ (3.5 - 1) / (6 - 1)
        @test Agrocosm.compute_lambda_npp(cd, 6.0) ≈ 1.0          # at LAI_max
        @test Agrocosm.compute_lambda_npp(cd, 9.0) == 1.0          # above LAI_max
    end

    @testset "litterfall and tendency (per second)" begin
        # Turnover is applied per second → tiny magnitudes.
        Λ = Agrocosm.compute_litterfall(cd, 1.0)
        @test Λ > 0
        @test Λ < 1.0e-6                                           # per-second, not per-year
        @test Agrocosm.compute_litterfall(cd, 2.0) ≈ 2 * Λ         # linear in LAI_b
        # dC/dt = (1-λ_NPP)·NPP − Λ_loc; below LAI_min λ_NPP=0 so growth ≈ NPP.
        NPP = 1.0e-8
        tend = Agrocosm.compute_carbon_tendency(cd, 0.5, NPP)      # LAI_b=0.5 < LAI_min
        @test tend ≈ NPP - Agrocosm.compute_litterfall(cd, 0.5)
    end
end
