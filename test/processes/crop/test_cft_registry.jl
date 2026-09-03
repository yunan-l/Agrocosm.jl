@testset "LPJmL crop CFT registry" begin
    @test length(CFTS) == 12
    @test getfield.(CFTS, :name) == Tuple(1:12)
    @test all(cft -> cft.plant_type == 2, CFTS)
    @test CFT_NAMES[4] == "tropical cereals"
    @test CFT_NAMES[9] == "oil crops soybean"
    @test crop_cft(4) === cft4
    @test crop_cft("oil crops soybean") === cft9

    @test cft1.path == 1
    @test cft1.laimax == 7.0f0
    @test cft1.hlimit == 360
    @test cft4.path == 2
    @test cft4.fphusen == 0.85f0
    @test cft9.path == 1
    @test cft9.basetemp.low == 7.0f0
    @test cft9.biological_fixation.enabled == 1
    @test cft9.biological_fixation.temperature_limit.low == 5.0f0
    @test cft9.biological_fixation.temperature_limit.high == 44.0f0
    @test cft1.biological_fixation.enabled == 0
    @test cft12.path == 2
    @test cft12.hiopt == 0.8f0

    expected_sowing_methods = Int32.(
        (2, 1, 4, 1, 4, 3, 1, 3, 1, 1, 2, 4),
    )
    expected_spring_temperatures = Float32.(
        (5, 18, 14, 12, 10, 8, 22, 13, 13, 15, 5, 14),
    )
    @test Tuple(cft.sowing_date.method for cft in CFTS) == expected_sowing_methods
    @test Tuple(cft.sowing_date.temp_spring for cft in CFTS) ==
          expected_spring_temperatures

    @test_throws ArgumentError crop_cft(13)
    @test_throws ArgumentError crop_cft("soybean")
end
