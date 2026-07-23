using Agrocosm
using Test

@testset "Crop soil mineralization / immobilization" begin
    m = CropNitrogenMineralization(Float64)   # soil_cn_ratio=15, immobilization_k=5e-3

    @testset "immobilization demand from C:N deficit" begin
        # litter_C/soil_CN − litter_N
        @test Agrocosm.immobilization_demand(m, 150.0, 5.0) ≈ 150.0 / 15 - 5.0   # = 5
        @test Agrocosm.immobilization_demand(m, 150.0, 20.0) == 0.0              # N-rich → net mineralization
    end

    @testset "Michaelis-Menten limitation" begin
        # concentration = available/depth·1000; limitation = conc/(k + conc) ∈ [0,1)
        lim = Agrocosm.immobilization_limitation(m, 1.0, 0.2)
        conc = 1.0 / 0.2 * 1000
        @test lim ≈ conc / (5.0e-3 + conc) rtol = 1e-12
        @test 0.0 ≤ lim < 1.0
        # more available nitrogen → less limiting
        @test Agrocosm.immobilization_limitation(m, 1.0, 0.2) > Agrocosm.immobilization_limitation(m, 1.0e-6, 0.2)
    end

    @testset "immobilized nitrogen capped at available" begin
        # ample supply → demand·limitation
        i1 = Agrocosm.immobilized_nitrogen(m, 5.0, 100.0, 0.2)
        @test i1 ≈ 5.0 * Agrocosm.immobilization_limitation(m, 100.0, 0.2) rtol = 1e-12
        @test i1 ≤ 5.0
        # scarce supply → capped at available
        @test Agrocosm.immobilized_nitrogen(m, 100.0, 0.5, 0.2) ≤ 0.5
        # no demand → no immobilization
        @test Agrocosm.immobilized_nitrogen(m, 0.0, 10.0, 0.2) == 0.0
    end
end
