using Agrocosm
using Test

@testset "Crop soil carbon dynamics" begin
    c = CropSoilCarbon(Float64)   # fast_fraction=0.98, atmospheric_fraction=0.5

    @testset "first-order decomposition" begin
        # (1 - exp(-rate·response))·pool
        @test Agrocosm.decomposed_carbon(0.1, 2.0, 100.0) ≈ (1 - exp(-0.2)) * 100.0
        @test Agrocosm.decomposed_carbon(0.1, 2.0, 0.0) == 0.0       # empty pool
        @test Agrocosm.decomposed_carbon(0.0, 2.0, 100.0) == 0.0     # no decomposition
        # decomposed fraction is in [0, pool] and monotone in rate·response
        @test Agrocosm.decomposed_carbon(1e6, 1.0, 100.0) ≈ 100.0 rtol = 1e-6   # fully decomposed
        @test Agrocosm.decomposed_carbon(0.2, 2.0, 100.0) > Agrocosm.decomposed_carbon(0.1, 2.0, 100.0)
    end

    @testset "litter routing conserves carbon" begin
        to_fast, to_slow, to_atm = Agrocosm.route_litter_carbon(c, 10.0)
        @test to_fast + to_slow + to_atm ≈ 10.0            # conservation
        @test to_atm ≈ 0.5 * 10.0                          # atmospheric fraction
        @test to_fast ≈ 0.98 * (0.5 * 10.0)                # 98% of the retained half
        @test to_slow ≈ 0.02 * (0.5 * 10.0)
        @test to_fast > to_slow
    end

    @testset "heterotrophic respiration" begin
        # atmospheric litter fraction + all decomposed fast/slow
        R = Agrocosm.heterotrophic_respiration(c, 10.0, 3.0, 1.0)
        @test R ≈ 0.5 * 10.0 + 3.0 + 1.0
    end
end
