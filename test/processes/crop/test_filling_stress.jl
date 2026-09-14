using Test
using Agrocosm
using Agrocosm: filling_weight, thermal_filling_progress, grain_sink_carbon,
                measured_filling_stress_exponent, cft3, cft1

@testset "filling is thermal until water is asked about" begin
    # The ablation contract: exponent 0 returns a weight of exactly one, so the
    # accumulated progress is the thermal progress bitwise.
    for sufficiency in (0.0, 0.31, 0.7, 1.0)
        @test filling_weight(sufficiency, 0.0) === 1.0
    end
    @test cft3.filling_stress_exponent == 0.0
    @test cft1.filling_stress_exponent == 0.0
    @test cft3.reserve_remobilisation == 1.0
end

@testset "a stressed day deposits less of its thermal time" begin
    # Monotone in stress and in the exponent, bounded by the unstressed day.
    @test filling_weight(1.0, 1.5) == 1.0
    @test filling_weight(0.0, 1.5) == 0.0
    @test filling_weight(0.6, 1.5) < filling_weight(0.8, 1.5) < 1.0
    @test filling_weight(0.6, 1.5) < filling_weight(0.6, 0.3)
    # Out-of-range sufficiency cannot make a grain heavier than its potential.
    @test filling_weight(1.4, 1.5) == 1.0
    @test filling_weight(-0.2, 1.5) == 0.0
end

@testset "the weighted integral is what sets grain weight" begin
    # Braunschweig maize: the dry arm's filling ran at a mean sufficiency that
    # the exponent turns into 0.837 of the potential grain weight. Integrating a
    # constant sufficiency over the whole window is the closed form of that.
    potential = grain_sink_carbon(3000.0, 0.126, 1.0, 1.0)
    @test potential ≈ 3000 * 0.126
    for (sufficiency, exponent) in ((0.6, 1.5), (0.8, 0.3), (0.9, 1.5))
        weighted = grain_sink_carbon(3000.0, 0.126, filling_weight(sufficiency, exponent), 1.0)
        @test weighted < potential
        @test weighted ≈ potential * sufficiency^exponent
    end
    # An unstressed season still fills every grain, which is what keeps the
    # mechanism silent where water is ample.
    @test grain_sink_carbon(3000.0, 0.126, filling_weight(1.0, 1.5), 1.0) ≈ potential
end

@testset "thermal progress is unchanged" begin
    @test thermal_filling_progress(0.70, 0.70) == 0.0
    @test thermal_filling_progress(1.00, 0.70) == 1.0
    @test thermal_filling_progress(0.85, 0.70) ≈ 0.5
    @test thermal_filling_progress(0.50, 0.70) == 0.0   # before the window closes
    @test thermal_filling_progress(0.80, 1.00) == 1.0   # a degenerate window
end

@testset "the measured exponents are per crop and per provenance" begin
    # Maize is measured against single-grain weight; wheat is fitted to a yield
    # contrast because no wheat deposit here recorded grain weight. Wheat must
    # come out markedly the less sensitive of the two - it remobilises 20-40% of
    # its grain carbon from stem reserves against maize's 10-20%.
    @test measured_filling_stress_exponent(3) == 5.0
    @test measured_filling_stress_exponent(1) == 3.0
    @test measured_filling_stress_exponent(3) > measured_filling_stress_exponent(1)
    # A crop no experiment has measured stays on thermal time.
    @test measured_filling_stress_exponent(2) == 0.0
    @test measured_filling_stress_exponent(9) == 0.0
end
