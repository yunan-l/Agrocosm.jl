using Agrocosm
using Test

@testset "Crop carbon allocation" begin
    a = CropCarbonAllocation(Float64)   # root_fraction_max=0.4, root_fraction_min=0.3

    @testset "root allocation fraction" begin
        # Early season (fphu = 0): maximum root fraction regardless of stress.
        @test Agrocosm.root_allocation_fraction(a, 0.0, 50.0) ≈ 0.4
        # Full stress (df = 0): stress term is 0 → still maximum root fraction.
        @test Agrocosm.root_allocation_fraction(a, 1.0, 0.0) ≈ 0.4
        # Well-watered late season: fewer roots (approaches max − min = 0.1 as df → large).
        f_late = Agrocosm.root_allocation_fraction(a, 1.0, 100.0)
        @test f_late < 0.4
        @test f_late ≈ 0.4 - 0.3 * (100.0 / (100.0 + exp(6.13 - 0.0883 * 100.0))) rtol = 1e-9
        # bounded in [max − min, max]
        for fphu in 0.0:0.1:1.0, df in 0.0:10.0:100.0
            fr = Agrocosm.root_allocation_fraction(a, fphu, df)
            @test 0.4 - 0.3 - 1e-9 ≤ fr ≤ 0.4 + 1e-9
        end
    end

    @testset "leaf carbon from LAI/SLA" begin
        # Constrained by LAI/SLA when carbon is ample.
        @test Agrocosm.leaf_carbon_from_lai(3.0, 0.03, 1000.0) ≈ 3.0 / 0.03
        # Constrained by available carbon when scarce.
        @test Agrocosm.leaf_carbon_from_lai(3.0, 0.03, 10.0) ≈ 10.0
        # Non-negative.
        @test Agrocosm.leaf_carbon_from_lai(3.0, 0.03, -5.0) == 0.0
    end
end
