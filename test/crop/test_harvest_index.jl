using Agrocosm
using Test

@testset "Crop harvest index" begin
    h = CropHarvestIndex(Float64)   # hiopt=0.5, himin=0.2
    hi(fphu, wdf) = Agrocosm.crop_harvest_index(h, fphu, wdf)

    fhiopt(fphu) = 100 * fphu / (100 * fphu + exp(11.1 - 10 * fphu))

    @testset "water-stress limits" begin
        # Well-watered (large wdf) → phenology-scaled optimum; fully stressed (wdf=0) → scaled minimum.
        @test hi(0.5, 1.0e6) ≈ fhiopt(0.5) * 0.5 rtol = 1e-6
        @test hi(0.5, 0.0) ≈ fhiopt(0.5) * 0.2 rtol = 1e-9
        # monotone increasing with water sufficiency
        vals = [hi(0.5, w) for w in 0.0:5.0:100.0]
        @test issorted(vals)
    end

    @testset "phenology dependence" begin
        # Grain filling: the optimal fraction rises with the heat-unit fraction.
        @test fhiopt(0.9) > fhiopt(0.5)
        @test hi(0.9, 1.0e6) > hi(0.5, 1.0e6)
    end

    @testset "closed form" begin
        f = fhiopt(0.5)
        hi_opt = f * 0.5
        hi_min = f * 0.2
        wf = 1.0 / (1.0 + exp(6.13 - 0.0883))
        @test hi(0.5, 1.0) ≈ (hi_opt - hi_min) * wf + hi_min rtol = 1e-9
    end

    @testset "HI above 1 scaled about 1" begin
        # e.g. sugarcane hiopt > 1: at full sufficiency, HI = fhiopt·(hiopt−1)+1.
        h2 = CropHarvestIndex(Float64; hiopt = 3.5, himin = 1.25)
        f = fhiopt(0.5)
        @test Agrocosm.crop_harvest_index(h2, 0.5, 1.0e6) ≈ f * (3.5 - 1) + 1 rtol = 1e-6
    end
end
