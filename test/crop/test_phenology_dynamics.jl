using Agrocosm
using Test

@testset "Crop phenological heat-unit accumulation" begin
    pd = CropPhenologyDynamics(Float64)   # heat_unit_requirement=1400, base_temperature=0

    @testset "heat-unit fraction" begin
        @test Agrocosm.heat_unit_fraction(pd, 0.0) == 0.0
        @test Agrocosm.heat_unit_fraction(pd, 700.0) ≈ 0.5
        @test Agrocosm.heat_unit_fraction(pd, 1400.0) ≈ 1.0
        @test Agrocosm.heat_unit_fraction(pd, 2000.0) == 1.0        # clamped at maturity
    end

    @testset "heat-unit accumulation rate (per second)" begin
        secs_per_day = 86400.0
        @test Agrocosm.heat_unit_rate(pd, 20.0) ≈ 20.0 / secs_per_day   # growing degree-days/day → per s
        @test Agrocosm.heat_unit_rate(pd, 0.0) == 0.0                   # at base temperature
        @test Agrocosm.heat_unit_rate(pd, -5.0) == 0.0                  # below base: no accumulation
        # integrating a constant 20 °C day for one day adds ~20 °C·day
        @test Agrocosm.heat_unit_rate(pd, 20.0) * secs_per_day ≈ 20.0
        # base temperature shifts the threshold
        pd8 = CropPhenologyDynamics(Float64; base_temperature = 8.0)
        @test Agrocosm.heat_unit_rate(pd8, 20.0) ≈ (20.0 - 8.0) / secs_per_day
        @test Agrocosm.heat_unit_rate(pd8, 5.0) == 0.0
    end
end
