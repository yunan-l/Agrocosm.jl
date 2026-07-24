using Agrocosm
using Test

@testset "Crop carbon pool budget" begin
    c = CropCarbon(Float64)   # specific_leaf_area = 30
    budget(biomass, fphu, β, lai, T, gpp) = Agrocosm.crop_carbon_budget(c, biomass, fphu, β, lai, T, gpp)

    @testset "organ partitioning conserves biomass" begin
        leaf, root, storage, npp = budget(1.0, 0.5, 1.0, 6.8, 25.0, 1.0e-7)
        @test leaf + root + storage ≈ 1.0            # organs sum to biomass
        @test all(≥(0), (leaf, root, storage))
        @test leaf ≈ 6.8 / 30 rtol = 1e-9            # leaf carbon set by LAI/SLA when carbon is ample
    end

    @testset "net primary production is GPP minus respiration" begin
        gpp = 1.0e-7
        _, _, _, npp = budget(1.0, 0.5, 1.0, 6.8, 25.0, gpp)
        @test 0.0 < npp < gpp                        # respiration reduces GPP but leaves a net gain
    end

    @testset "no biomass, no organs" begin
        leaf, root, storage, npp = budget(0.0, 0.5, 1.0, 6.8, 25.0, 1.0e-7)
        @test (leaf, root, storage) == (0.0, 0.0, 0.0)
        # no biomass → no maintenance respiration, but growth respiration still takes r_growth of GPP
        @test npp ≈ (1.0 - 0.25) * 1.0e-7 rtol = 1e-9
    end

    @testset "water stress shifts carbon to roots" begin
        _, root_wet, _, _ = budget(1.0, 0.5, 1.0, 6.8, 25.0, 1.0e-7)    # β = 1 (well watered)
        _, root_dry, _, _ = budget(1.0, 0.5, 0.0, 6.8, 25.0, 1.0e-7)    # β = 0 (full stress)
        @test root_dry ≥ root_wet
    end
end
