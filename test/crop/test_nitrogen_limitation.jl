using Agrocosm
using Test

@testset "Crop nitrogen limitation of Vcmax" begin
    lim = CropNitrogenVcmaxLimit(Float64)   # ncleaf_min=1/58.8, k_temp=0.0693, pressure_scale=25

    nlv(pot, N, C, T) = Agrocosm.nitrogen_limited_vcmax(lim, pot, N, C, T)

    @testset "non-limiting: abundant nitrogen" begin
        # Plenty of leaf N above the structural minimum → capacity exceeds potential → no limitation.
        v, f = nlv(10.0, 5.0, 100.0, 25.0)
        @test v ≈ 10.0
        @test f ≈ 1.0
    end

    @testset "limiting: scarce nitrogen caps Vcmax" begin
        # Just above structural minimum (1/58.8·100 ≈ 1.7007) → small Rubisco N → capacity < potential.
        v, f = nlv(10.0, 1.72, 100.0, 25.0)
        @test v < 10.0
        @test 0.0 < f < 1.0
        # matches the closed form: capacity = rubisco_N/(25e-3)·(86400·12·1e-6) at T=25
        rubisco_N = 1.72 - (1 / 58.8) * 100.0
        capacity = rubisco_N / (25.0e-3) * (86400 * 12 * 1.0e-6)
        @test v ≈ capacity rtol = 1e-9
        @test f ≈ capacity / 10.0 rtol = 1e-9
    end

    @testset "structural nitrogen protected" begin
        # Below the structural minimum → no Rubisco nitrogen → essentially zero capacity.
        v, f = nlv(10.0, 1.0, 100.0, 25.0)   # 1.0 < 1.7007 structural
        @test v ≤ 1e-6
        @test f ≈ 0.0 atol = 1e-6
    end

    @testset "no potential capacity" begin
        v, f = nlv(0.0, 5.0, 100.0, 25.0)
        @test v == 0.0
        @test f == 0.0
    end

    @testset "temperature dependence" begin
        # Higher temperature raises the nitrogen requirement (exp(-k(T-25)) < 1 in the denominator
        # increases capacity), so warmer → higher N-limited capacity for the same leaf N.
        _, f_cold = nlv(10.0, 1.72, 100.0, 15.0)
        _, f_warm = nlv(10.0, 1.72, 100.0, 35.0)
        @test f_warm > f_cold
    end
end
