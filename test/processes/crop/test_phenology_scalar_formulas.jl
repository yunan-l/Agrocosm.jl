using Agrocosm
using Test

@testset "Phenology scalar formulas" begin
    T = Float32
    effective = (low = T(0), high = T(15))
    optimum = (low = T(5), high = T(10))

    @test Agrocosm.compute_phenology_fraction(T(20), T(40)) == T(0.5)
    @test Agrocosm.compute_phenology_fraction(T(20), T(0)) == T(0)
    @test Agrocosm.compute_heat_unit_increment(T(3), T(5)) == T(0)
    @test Agrocosm.compute_heat_unit_increment(T(8), T(5)) == T(3)

    @test Agrocosm.compute_vernalization_increment(T(2.5), T(0), T(20), effective, optimum) == T(0.5)
    @test Agrocosm.compute_vernalization_increment(T(7), T(0), T(20), effective, optimum) == T(1)
    @test Agrocosm.compute_vernalization_increment(T(12.5), T(0), T(20), effective, optimum) == T(0.5)
    @test Agrocosm.compute_vernalization_increment(T(7), T(20), T(20), effective, optimum) == T(0)

    @test Agrocosm.compute_vernalization_factor(T(3), T(20)) == T(0)
    @test Agrocosm.compute_vernalization_factor(T(12), T(20)) == T(0.5)
    @test Agrocosm.compute_vernalization_factor(T(20), T(20)) == T(1)

    @test Agrocosm.compute_photoperiod_factor(T(0.2), T(0.7), T(0.2), T(8), T(10), T(14)) == T(0.2)
    @test Agrocosm.compute_photoperiod_factor(T(0.2), T(0.7), T(0.2), T(14), T(10), T(14)) == T(1)
    @test Agrocosm.compute_photoperiod_factor(T(0.8), T(0.7), T(0.2), T(8), T(10), T(14)) == T(1)

    presenescent = Agrocosm.compute_phenology_lai_fraction(
        T(0.4), T(0.05), T(0.05), T(0.45), T(0.45), T(0.7), T(0), T(1),
    )
    senescent = Agrocosm.compute_phenology_lai_fraction(
        T(0.8), T(0.05), T(0.05), T(0.45), T(0.45), T(0.7), T(0), T(1),
    )
    mature = Agrocosm.compute_phenology_lai_fraction(
        T(1), T(0.05), T(0.05), T(0.45), T(0.45), T(0.7), T(0), T(0.5),
    )
    @test zero(T) < presenescent < one(T)
    @test senescent ≈ T(2 / 3) atol = eps(T)
    @test mature == zero(T)
end
