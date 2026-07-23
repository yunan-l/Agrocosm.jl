using Agrocosm
using Test

@testset "Crop nitrogen allocation" begin
    a = CropNitrogenAllocation(Float64)   # ratio_root=1.16, ratio_storage=0.99, ratio_pool=3
    alloc(N, lc, rc, sc, pc) = Agrocosm.allocate_crop_nitrogen(a, N, lc, rc, sc, pc)

    @testset "conserves total nitrogen" begin
        ln, rn, sn, pn = alloc(10.0, 100.0, 50.0, 5.0, 10.0)
        @test ln + rn + sn + pn ≈ 10.0
        @test all(≥(0), (ln, rn, sn, pn))
    end

    @testset "weights by carbon / target C:N ratio" begin
        ln, rn, sn, pn = alloc(10.0, 100.0, 50.0, 5.0, 10.0)
        # leaf:root nitrogen ratio = leaf_c / (root_c/ratio_root)
        @test ln / rn ≈ 100.0 / (50.0 / 1.16) rtol = 1e-9
        @test ln / pn ≈ 100.0 / (10.0 / 3.0) rtol = 1e-9
    end

    @testset "leaf-only crop" begin
        ln, rn, sn, pn = alloc(4.0, 100.0, 0.0, 0.0, 0.0)
        @test ln ≈ 4.0
        @test (rn, sn, pn) == (0.0, 0.0, 0.0)
    end

    @testset "no organ carbon" begin
        @test alloc(10.0, 0.0, 0.0, 0.0, 0.0) == (0.0, 0.0, 0.0, 0.0)
    end
end
