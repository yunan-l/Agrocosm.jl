using Test
using Agrocosm
using Agrocosm: establishment_loss, measured_establishment_loss_rate, cft1, cft3

@testset "nothing is lost until a rate is asked for" begin
    # The ablation contract: at rate 0 no day destroys any stand, however hot.
    for temperature in (20.0, 35.0, 45.0, 50.0)
        @test establishment_loss(temperature, 35.0, 0.0, 0.3) === 0.0
    end
    @test cft1.establishment_loss_rate == 0.0
    @test cft3.establishment_loss_rate == 0.0
end

@testset "the ceiling is wheat's published germination maximum" begin
    # Porter and Gawith (1999) give wheat germination a maximum of 35 C. The
    # loss is negligible below it and near the full rate well above.
    @test cft1.establishment_heat_ceiling == 35.0
    @test establishment_loss(30.0, 35.0, 0.06, 0.3) < 0.06 * 1e-6
    @test establishment_loss(35.0, 35.0, 0.06, 0.3) ≈ 0.06 / 2
    @test establishment_loss(40.0, 35.0, 0.06, 0.3) ≈ 0.06 rtol = 1e-6
    # Monotone in temperature and linear in the rate.
    @test establishment_loss(36.0, 35.0, 0.06, 0.3) > establishment_loss(35.5, 35.0, 0.06, 0.3)
    @test establishment_loss(40.0, 35.0, 0.12, 0.3) ≈ 2 * establishment_loss(40.0, 35.0, 0.06, 0.3)
end

@testset "the window is days, not thermal time" begin
    # The reason this matters is measurable: a window closed on `fphu` closes
    # fastest in the hottest season. Across the Hot Serial Cereal sowings the
    # thermal window held open 5 days for the three that failed in the field and
    # 13 for the January sowings that yielded 8-10 t/ha, so the crop escaped the
    # stress the window exists to represent.
    @test cft1.establishment_days == 30
    @test cft1.establishment_days isa Integer
end

@testset "the rate was measured against failures, on both sides" begin
    # 0.06 is where three field failures reach exactly zero AND the twelve
    # harvested treatments do not move. A crop no experiment has failed gets
    # nothing, which leaves it unable to fail - as the whole lineage is.
    @test measured_establishment_loss_rate(1) == 0.06
    @test measured_establishment_loss_rate(3) == 0.0
    @test measured_establishment_loss_rate(2) == 0.0
end

@testset "flowering heat is what reconciles the two sites" begin
    # `heat_day_rate` ships at 0.02, which no experiment had constrained. Hot
    # Serial Cereal flowered at 39 C and says 0.008: at that rate the Maricopa
    # grain-number bias is 1.05 against 1.48 with the mechanism off, while
    # Braunschweig - which never reaches the 38 C threshold - is untouched.
    @test Agrocosm.measured_anthesis_heat_rate(1) == 0.008
    @test Agrocosm.measured_anthesis_heat_rate(3) == 0.0
    # Wheat overrides the generic 38 C threshold with its own 36.
    @test Agrocosm.cft1.heat_day_temperature == 36.0
    # The shipped rate is the one this measurement corrects, so the two must
    # differ and the measurement must be the smaller.
    @test Agrocosm.measured_anthesis_heat_rate(1) < Agrocosm.cft1.heat_day_rate
end
