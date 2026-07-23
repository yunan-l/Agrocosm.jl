using Agrocosm
using Terrarium
using Test

@testset "Crop root distribution (LPJmL exponential)" begin
    β = 0.94
    rd = CropRootDistribution(Float64; beta_root = β)

    @testset "continuous density" begin
        # ∂R/∂z ∝ β^{d_cm}, d_cm = -100 z; carries the -ln β prefactor.
        @test Terrarium.root_density(rd, 0.0) ≈ -log(β)
        @test Terrarium.root_density(rd, -0.10) ≈ -log(β) * β^10   # 10 cm depth
        @test Terrarium.root_density(rd, -0.50) ≈ -log(β) * β^50   # 50 cm depth
        # density decreases with depth (more negative z)
        @test Terrarium.root_density(rd, -0.5) < Terrarium.root_density(rd, -0.1)
    end

    @testset "normalized column fractions" begin
        # Fine uniform grid, 2 cm layers to 3 m depth (β^300 ≈ 0 so the column
        # captures essentially all roots).
        grid = ColumnGrid(CPU(), Float64, UniformSpacing(Δz = 0.02, N = 150))
        rf = Field(Agrocosm.crop_root_fraction(rd, grid, nothing, (;)))
        compute!(rf)
        frac = interior(rf)[1, 1, :]                 # bottom → top ordering
        z = znodes(get_field_grid(grid), Center(), Center(), Center())

        @test all(>(0), frac)
        @test sum(frac) ≈ 1.0                        # normalized
        # roots concentrate near the surface: fraction increases toward z = 0
        @test issorted(frac)                         # ascending with z (surface last)

        # Cumulative fraction above a depth approaches the analytic LPJmL CDF
        # Y(d) = 1 - β^{d_cm}. Compare at ~50 cm depth.
        target_depth_cm = 50.0
        above = sum(frac[i] for i in eachindex(frac) if -100 * z[i] ≤ target_depth_cm)
        @test above ≈ (1 - β^target_depth_cm) atol = 0.02
    end
end
