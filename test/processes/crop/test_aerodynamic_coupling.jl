using Test
using Agrocosm
using Agrocosm: aerodynamic_conductance, coupled_demand_parameters,
                compute_transpiration_demand, saturation_vapour_pressure,
                saturation_vapour_pressure_slope, fao56_canopy_height, lpjmlparams

@testset "the coupling is off until it is asked for" begin
    # The ablation contract: at rate 0 the two constants come back untouched, so
    # every earlier run reproduces bitwise.
    alpha, shape = coupled_demand_parameters(
        1.485, 2.2, 3.0, 22.0, 1000.0, 3.46, 2.0, 15.0, 0.0,
    )
    @test alpha === 1.485
    @test shape === 2.2
    # And with the rate on but no humidity in the forcing, which is how a run
    # without `huss` must behave rather than inventing a deficit.
    @test coupled_demand_parameters(1.485, 2.2, 3.0, 22.0, 0.0, 3.46, 2.0, 15.0, 1.0) ===
          (1.485, 2.2)
    # A bare-soil day, where the height scaling has driven the canopy to zero.
    @test coupled_demand_parameters(1.485, 2.2, 3.0, 22.0, 1000.0, 3.46, 0.0, 15.0, 1.0) ===
          (1.485, 2.2)
end

@testset "the equivalent parameters reproduce Penman-Monteith" begin
    # The whole point of the substitution: LPJmL's demand evaluated with these
    # two numbers must equal PM evaluated independently, to round-off. If this
    # fails the substitution is not an identity and the mechanism is a fit.
    for (eeq, temp, deficit, wind, height, gc, daylength) in
        ((3.0, 22.0, 1000.0, 3.46, 2.0, 10.0, 15.0), (2.0, 18.0, 500.0, 1.7, 0.9, 6.0, 11.0),
         (5.0, 30.0, 2500.0, 1.7, 0.9, 14.0, 13.5), (1.5, 12.0, 300.0, 5.0, 1.6, 4.0, 9.0))
        alpha, shape = coupled_demand_parameters(
            1.485, 2.2, eeq, temp, deficit, wind, height, daylength, 1.0,
        )
        through_lpjml = compute_transpiration_demand(0.0, eeq, alpha, shape, gc)
        ga = aerodynamic_conductance(wind, height, 2.0) / 1000
        ratio = saturation_vapour_pressure_slope(temp) / 66.5 + 1
        # The imposed term is daylight-weighted to match `eeq`'s own convention.
        penman = (eeq * ratio + 0.6477 * deficit * ga * daylength / 24) /
                 (ratio + ga / (gc / 1000))
        @test through_lpjml ≈ penman rtol = 1e-10
    end
end

@testset "coupling is a canopy property, not a constant" begin
    # LPJmL's sensitivity to stomatal closure is the same number everywhere
    # because GM is fixed; the measured one is not. Maize in wind must come out
    # markedly more sensitive than wheat in still air.
    sensitivity(alpha, shape, gc) = shape * alpha / (gc + shape * alpha)
    maize = coupled_demand_parameters(1.485, 2.2, 3.0, 22.0, 1000.0, 3.46, 2.0, 15.0, 1.0)
    wheat = coupled_demand_parameters(1.485, 2.2, 3.0, 22.0, 1000.0, 1.7, 0.9, 15.0, 1.0)
    @test sensitivity(maize..., 10.0) > sensitivity(wheat..., 10.0)
    @test sensitivity(1.485, 2.2, 10.0) < sensitivity(wheat..., 10.0)
    # Taller and windier both tighten the coupling.
    @test aerodynamic_conductance(3.46, 2.0, 2.0) > aerodynamic_conductance(1.7, 2.0, 2.0)
    @test aerodynamic_conductance(1.7, 2.0, 2.0) > aerodynamic_conductance(1.7, 0.9, 2.0)
end

@testset "the deficit must not come from the daily mean temperature" begin
    # FAO-56 requires the saturation pressure as the mean of es(Tmax) and
    # es(Tmin). Using es(Tmean) collapses the deficit on humid days, and in this
    # model that took Braunschweig's season transpiration from 237 mm to 56.
    mean_temperature, range = 17.0, 12.0
    deficits = map((1500.0, 1800.0)) do actual
        from_mean = saturation_vapour_pressure(mean_temperature) - actual
        fao56 = (saturation_vapour_pressure(mean_temperature + range / 2) +
                 saturation_vapour_pressure(mean_temperature - range / 2)) / 2 - actual
        (from_mean, fao56)
    end
    @test all(d -> d[2] > d[1], deficits)
    # The gap widens as the air approaches saturation, which is exactly the
    # regime a humid site spends its season in: 1.28 at a middling day, 1.9 on a
    # near-saturated one.
    @test deficits[1][2] / deficits[1][1] > 1.25
    @test deficits[2][2] / deficits[2][1] > 1.8
end

@testset "canopy heights are the published ones" begin
    @test fao56_canopy_height(3) == 2.0    # maize, FAO-56 Table 12
    @test fao56_canopy_height(1) == 1.0    # wheat
    @test fao56_canopy_height(9) == 0.75   # soybean
    # A crop the table does not cover stays on the uncoupled demand.
    @test fao56_canopy_height(7) == 0.0
    @test lpjmlparams.aerodynamic_coupling == 0.0
end
