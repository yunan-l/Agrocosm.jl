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

@testset "measured senescence shape" begin
    # Maricopa FACE logged leaf area through senescence on four irrigation arms.
    # Solving k in ((1 - fphu) / (1 - fphusen))^k from each measurement gives
    # 2.10 and 1.67 on the two DRY arms and 0.31 and 0.51 on the two WET ones, so
    # the shipped 2 is the water-stressed crop's value. Only wheat was measured.
    @test Agrocosm.measured_senescence_shape(1) == 0.41
    @test Agrocosm.measured_senescence_shape(3) == 0.0
    @test Agrocosm.cft1.shapesenescencenorm == 2

    # What the difference is worth: at fphu 0.89, two thirds of the way through
    # senescence, the shipped curve has taken the canopy to an eighth of its peak
    # and the measured one holds two thirds. That gap is the model's missing
    # transpiration in the hottest weeks of the season.
    shipped = Agrocosm.compute_phenology_lai_fraction(
        0.89, 0.05, 0.05, 0.45, 0.45, 0.70, 0.0, 2.0,
    )
    measured = Agrocosm.compute_phenology_lai_fraction(
        0.89, 0.05, 0.05, 0.45, 0.45, 0.70, 0.0,
        Agrocosm.measured_senescence_shape(1),
    )
    @test shipped ≈ 0.1344 atol = 1e-4
    @test measured ≈ 0.6628 atol = 1e-4
    @test measured / shipped > 4
end
