@testset "LPJmL crop PFT registry" begin
    @test length(CROP_PFTS) == 12
    @test getfield.(CROP_PFTS, :name) == Tuple(1:12)
    @test CROP_PFT_NAMES[4] == "tropical cereals"
    @test CROP_PFT_NAMES[9] == "oil crops soybean"
    @test crop_pft(4) === cft4
    @test crop_pft("oil crops soybean") === cft9

    @test cft1.path == 1
    @test cft1.laimax == 7.0f0
    @test cft1.hlimit == 360
    @test cft4.path == 2
    @test cft4.fphusen == 0.85f0
    @test cft9.path == 1
    @test cft9.basetemp.low == 7.0f0
    @test cft12.path == 2
    @test cft12.hiopt == 0.8f0

    @test_throws ArgumentError crop_pft(13)
    @test_throws ArgumentError crop_pft("soybean")
end
