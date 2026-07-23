using Agrocosm
using Test

@testset "Crop phenology (LPJmL LAI trajectory)" begin
    p = CropPhenology(Float64)   # fphuc=.05,flaimaxc=.05,fphuk=.45,flaimaxk=.95,fphusen=.70,laimax=7

    @testset "trajectory anchors" begin
        # By construction the growth curve passes through (fphuc, flaimaxc) and (fphuk, flaimaxk),
        # and the senescence branch starts at 1 at fphusen and reaches 0 at fphu = 1.
        @test Agrocosm.compute_lai_fraction(p, 0.0) ≈ 0.0 atol = 1e-12
        @test Agrocosm.compute_lai_fraction(p, 0.05) ≈ 0.05 rtol = 1e-6   # flaimaxc at fphuc
        @test Agrocosm.compute_lai_fraction(p, 0.45) ≈ 0.95 rtol = 1e-6   # flaimaxk at fphuk
        @test Agrocosm.compute_lai_fraction(p, 0.70) ≈ 1.0 rtol = 1e-6    # peak at senescence onset
        @test Agrocosm.compute_lai_fraction(p, 1.0) ≈ 0.0 atol = 1e-12    # bare at end of cycle
        # senescence: ((1-fphu)/(1-fphusen))^shape, shape=2 → at 0.85, (0.15/0.30)^2 = 0.25
        @test Agrocosm.compute_lai_fraction(p, 0.85) ≈ 0.25 rtol = 1e-6
    end

    @testset "LAI scaling and shape" begin
        @test Agrocosm.compute_crop_lai(p, 0.70) ≈ 7.0 rtol = 1e-6        # peak LAI = laimax
        @test Agrocosm.compute_crop_lai(p, 0.0) ≈ 0.0 atol = 1e-12
        # rises to the senescence onset, then declines
        rising = [Agrocosm.compute_lai_fraction(p, f) for f in 0.0:0.05:0.70]
        @test issorted(rising)
        falling = [Agrocosm.compute_lai_fraction(p, f) for f in 0.70:0.05:1.0]
        @test issorted(falling; rev = true)
        # bounded in [0, 1] across the whole cycle
        for f in 0.0:0.02:1.0
            frac = Agrocosm.compute_lai_fraction(p, f)
            @test 0.0 ≤ frac ≤ 1.0 + 1e-9
        end
    end

    @testset "harvest floor" begin
        # With a non-zero harvest floor, LAI does not decline all the way to zero.
        ph = CropPhenology(Float64; flaimaxharvest = 0.5)
        @test Agrocosm.compute_lai_fraction(ph, 1.0) ≈ 0.5 rtol = 1e-6
    end
end
