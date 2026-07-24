using Agrocosm
using Test
using Enzyme
using FiniteDifferences

# Differentiability of the crop scalar primitives. Terrarium targets full (Enzyme) differentiability;
# these check that the ported crop physics functions carry valid reverse-mode adjoints on the CPU —
# each gradient is compared against a finite-difference reference. (The whole-model rollout is
# differentiated through Reactant on a GPU-capable build; see docs/dev/2026-07/spike_crop_reactant_ad.jl.)
@testset "Crop primitive differentiability (Enzyme reverse mode)" begin
    fdm = central_fdm(5, 1)
    reverse_grad(f, x) = Enzyme.autodiff(Reverse, f, Active, Active(x))[1][1]

    @testset "net primary production wrt gross assimilation" begin
        growth = CropGrowthRespiration(Float64)
        # Well inside the positive (growth) regime so the finite-difference stencil stays clear of the
        # kink at gross = maintenance; the adjoint is the constant growth-respiration factor (1−r_growth).
        f(gpp) = net_primary_production(growth, gpp, 0.1)
        @test reverse_grad(f, 1.0) ≈ fdm(f, 1.0) rtol = 1e-6
        @test 0 < reverse_grad(f, 1.0) < 1
    end

    @testset "soil decomposition temperature response wrt temperature" begin
        response = CropSoilDecompositionResponse(Float64)
        # The combined response is clamped to [0,1] (flat where it saturates); differentiate the
        # unclamped Lloyd-Taylor temperature response, which is smooth and increasing.
        f(T) = soil_decomposition_temperature_response(response, T)
        @test reverse_grad(f, 5.0) ≈ fdm(f, 5.0) rtol = 1e-6
        @test reverse_grad(f, 5.0) > 0   # warmer soil decomposes faster
    end

    @testset "heterotrophic respiration wrt litter carbon" begin
        bgc = CropSoilBiogeochemistry(Float64)
        f(litter) = soil_carbon_tendencies(bgc, litter, 5.0, 20.0, 1.0)[4]   # het component
        @test reverse_grad(f, 1.0) ≈ fdm(f, 1.0) rtol = 1e-6
    end

    @testset "gross nitrification wrt ammonium" begin
        nitrification = CropNitrification(Float64)
        f(ammonium) = gross_nitrification(nitrification, ammonium, 0.6, 25.0, 6.5)[1]
        @test reverse_grad(f, 0.1) ≈ fdm(f, 0.1) rtol = 1e-6
    end

    @testset "leaf nitrogen limitation wrt leaf nitrogen (linear regime)" begin
        nitrogen = CropNitrogen(Float64)
        # Choose a leaf N:C between ncleaf_min and ncleaf_ref so the limitation is in its linear regime.
        leaf_carbon = 1.0
        f(leaf_nitrogen) = leaf_nitrogen_limitation(nitrogen, leaf_nitrogen, leaf_carbon)
        leaf_nitrogen = 0.5 * (nitrogen.ncleaf_min + nitrogen.ncleaf_ref) * leaf_carbon
        @test reverse_grad(f, leaf_nitrogen) ≈ fdm(f, leaf_nitrogen) rtol = 1e-6
        @test reverse_grad(f, leaf_nitrogen) > 0   # more leaf N relaxes the limitation
    end
end
