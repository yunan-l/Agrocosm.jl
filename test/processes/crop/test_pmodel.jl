using Test
using Agrocosm

# The P-model, validated against the published reference implementation.
#
# This is the one mechanism in this project with an external ground truth: the
# rpmodel vignette publishes intermediate values for a stated case, so the physics
# can be checked before any of it is wired into the model. Everything else here
# has had to be argued from first principles or measured against yields.

const TM = Float64

@testset "water viscosity matches reference values" begin
    # The least-cost ratio carries `1.6 * eta_star`, so a wrong viscosity moves
    # the optimal ci/ca and nothing else complains. It is checked against water
    # rather than against the P-model for exactly that reason: laid out
    # density-major, Huber's double sum was indexed the wrong way round here and
    # gave 1e7 Pa s, which shifted xi by only 4% and would have been easy to miss.
    for (celsius, reference) in ((0.0, 1.7914e-3), (10.0, 1.3060e-3),
                                 (20.0, 1.0016e-3), (25.0, 0.8900e-3),
                                 (30.0, 0.7972e-3), (40.0, 0.6527e-3))
        @test Agrocosm.pmodel_viscosity_ratio(TM(celsius), TM(101325)) ≈
              TM(reference) rtol = 1e-3
    end
end

@testset "the published rpmodel case is reproduced" begin
    # rpmodel vignette: 20 C, sea level, VPD 1000 Pa, 400 ppm.
    tc, patm, vpd = TM(20), TM(101325), TM(1000)
    ca = TM(400) * TM(1e-6) * patm
    optimum = Agrocosm.pmodel_optimal_chi(tc, patm, vpd, ca, TM(146))
    @test optimum.kmm ≈ 46.09928 rtol = 1e-6
    @test optimum.gammastar ≈ 3.339251 rtol = 1e-6
    @test optimum.xi ≈ 63.3145 rtol = 1e-5
    @test optimum.chi ≈ 0.694352 rtol = 1e-5
end

@testset "the optimality responses go the right way" begin
    patm, ca = TM(101325), TM(400) * TM(1e-6) * TM(101325)
    chi(tc, vpd) = Agrocosm.pmodel_optimal_chi(TM(tc), patm, TM(vpd), ca, TM(146)).chi
    # A thirstier atmosphere closes stomata: ci/ca falls with vapour deficit, and
    # as a square root rather than linearly - the signature of the theory.
    @test chi(20, 2000) < chi(20, 1000) < chi(20, 500)
    # Warmer air RAISES the optimal ci/ca, which is the opposite of the naive
    # expectation and is the interesting half of the theory: the Michaelis-Menten
    # coefficient climbs Arrhenius-fast with temperature while water viscosity
    # falls, so `xi` rises and stomata can afford to stay open. Measured here:
    # chi runs 0.486 at 5 C to 0.885 at 40 C. A prescribed lambda cannot express
    # that, which is one of the things this path is being tested for.
    @test chi(30, 1000) > chi(10, 1000)
    @test chi(5, 1000) < chi(20, 1000) < chi(40, 1000)
    # Zero deficit is the no-cost limit, not a division by zero.
    @test 0 < chi(20, 0) <= 1
end

@testset "Jmax limitation is bounded and guarded at its singular point" begin
    gammastar = TM(3.339251)
    ca = TM(400) * TM(1e-6) * TM(101325)
    m(chi) = Agrocosm.pmodel_limitation(TM(chi), ca, gammastar, TM(0.41))
    @test 0 < m(0.7) < 1
    @test m(0.9) > m(0.5)
    # Below `cstar` the root turns imaginary; the guard returns zero rather than
    # NaN, which matters because this file is on the differentiated path.
    @test m(0.05) == zero(TM)
    @test isfinite(m(0.0))
end

@testset "C4 drops the photorespiration term rather than special-casing it" begin
    tc, patm, vpd = TM(25), TM(101325), TM(1500)
    ca = TM(400) * TM(1e-6) * patm
    c3 = Agrocosm.pmodel_light_use_efficiency(tc, patm, vpd, ca, TM(0.05), TM(146), TM(0.41))
    c4 = Agrocosm.pmodel_light_use_efficiency(tc, patm, vpd, ca, TM(0.05), TM(146), TM(0))
    # C4 concentrates CO2 at Rubisco, so `m'` is one and efficiency is higher.
    @test c4 > c3 > 0
    @test c4 ≈ TM(0.05) * Agrocosm.pmodel_quantum_yield_factor(tc) * TM(12.0107)
end

@testset "vapour pressure deficit from specific humidity and pressure" begin
    vpd(tc, patm, q) = Agrocosm.pmodel_vapour_pressure_deficit(TM(tc), TM(patm), TM(q))
    # Saturated air has no deficit; the conversion must not go negative there.
    saturated_q(tc, patm) = begin
        es = 611.2 * exp(17.62 * tc / (243.12 + tc))
        0.622 * es / (patm - 0.378 * es)
    end
    @test vpd(20, 101325, saturated_q(20, 101325)) ≈ 0 atol = 1e-6
    @test vpd(20, 101325, 0.0) ≈ 611.2 * exp(17.62 * 20 / (243.12 + 20)) rtol = 1e-9
    # Drier air, larger deficit; warmer air at the same humidity, larger deficit.
    @test vpd(20, 101325, 0.005) > vpd(20, 101325, 0.010)
    @test vpd(30, 101325, 0.008) > vpd(20, 101325, 0.008)
    # Negative humidity is data error, not a negative deficit.
    @test vpd(20, 101325, -0.001) > 0
end
