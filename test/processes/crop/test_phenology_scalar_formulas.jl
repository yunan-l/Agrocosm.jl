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

@testset "two-exponential root profile" begin
    # Zero on either rate is the single exponential, bitwise: the ablation.
    @test Agrocosm.root_distribution(0.94) == Agrocosm.root_distribution(0.94, 0.0, 0.0)
    @test Agrocosm.root_distribution(0.94) == Agrocosm.root_distribution(0.94, 0.13, 0.0)

    surface, deep = Agrocosm.measured_root_profile(1)
    measured = Agrocosm.root_distribution(0.94, surface, deep)
    shipped = Agrocosm.root_distribution(0.94)
    @test sum(measured) ≈ 1 atol = 1e-6
    # Maricopa cored 0.148 of root mass below 45 cm; the single exponential
    # carries 0.062 there and puts 0.21% below a metre.
    @test sum(shipped[4:5]) ≈ 0.0021 atol = 1e-4
    @test sum(measured[4:5]) ≈ 0.0335 atol = 1e-3
    @test sum(measured[4:5]) / sum(shipped[4:5]) > 15
    # It does not achieve that by emptying the surface: the cores put 0.598 in
    # 0-15 cm, so the top layer must stay dominant.
    @test measured[1] > 0.65

    @test Agrocosm.measured_root_profile(3) == (0.0, 0.0)
end

@testset "expansion water stress" begin
    # Recorded, not installed: the rate is zero until the reach question that
    # Braunschweig 2015 raises is resolved, so it cannot move a shipped result.
    threshold, rate = Agrocosm.measured_expansion_stress(1)
    @test rate == 0.0
    # A matric potential in bar, inside wheat's published -0.5 to -1.0 onset.
    @test 0.5 < threshold < 1.0

    # The weight is one when the crop is wetter than the threshold, and falls
    # log-linearly to zero at the permanent wilting point.
    @test Agrocosm.expansive_growth_weight(-0.40, threshold) == 1.0
    @test Agrocosm.expansive_growth_weight(-15.0, threshold) == 0.0
    # Maricopa's deficit arm: -1.22 bar, and its canopy thermometer measured
    # that arm transpiring at 0.80 of its fully irrigated one.
    @test Agrocosm.expansive_growth_weight(-1.22, threshold) ≈ 0.80 atol = 0.01
    # Threshold zero returns one, which is the ablation.
    @test Agrocosm.expansive_growth_weight(-5.0, 0.0) === 1.0

    # Campbell through the model's own two anchors, -15 bar and -1/3 bar.
    @test Agrocosm.soil_water_potential(1.0, 0.197, 0.316) ≈ -1 / 3 atol = 1e-6
    @test Agrocosm.soil_water_potential(0.0, 0.197, 0.316) ≈ -15.0 atol = 0.1
    # The whole point: the SAME relative content is a different potential in a
    # clay loam and a loamy sand, and the crop's thermometer ranks them the way
    # the potential does.
    clay = Agrocosm.soil_water_potential(0.6, 0.197, 0.316)
    sand = Agrocosm.soil_water_potential(0.6, 0.048, 0.214)
    @test clay < sand

    @test Agrocosm.expansion_stress_loss(0.8, 1.0, true, 0.05) ≈ 0.01 atol = 1e-9
    @test Agrocosm.expansion_stress_loss(1.0, 1.0, true, 0.05) == 0.0
    @test Agrocosm.expansion_stress_loss(0.8, 1.0, false, 0.05) == 0.0
end

@testset "grain filling temperature weight" begin
    optimum, rate = Agrocosm.measured_filling_temperature(1)
    # Rate zero returns exactly one, so the product with the thermal increment is
    # bitwise the inherited filling. That is the ablation contract.
    @test Agrocosm.filling_temperature_weight(40.0, optimum, 0.0) === 1.0
    @test Agrocosm.filling_temperature_weight(-10.0, optimum, 0.0) === 1.0
    # Below the optimum the grain still cannot exceed its shipped potential: the
    # weight only ever subtracts, like every other mechanism here.
    @test Agrocosm.filling_temperature_weight(optimum - 10, optimum, rate) == 1.0
    @test Agrocosm.filling_temperature_weight(optimum, optimum, rate) == 1.0
    # And it is bounded below, so an extreme day cannot make the sink negative.
    @test Agrocosm.filling_temperature_weight(200.0, optimum, rate) == 0.0
    # The measured slope: 0.0855 of a day's filling per degree above 27 C.
    @test Agrocosm.filling_temperature_weight(optimum + 5, optimum, rate) ≈
          1 - 5 * rate atol = 1e-12
    @test Agrocosm.measured_filling_temperature(3) == (0.0, 0.0)

    # The shipped potential is 40.0 mg of dry matter at 0.45 gC/g.
    @test Agrocosm.published_grain_traits(1)[2] * 1000 / 0.45 ≈ 40.0 atol = 0.1
end
