using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA prescribed fertilizer is conserved" begin
    crop, crop_cal, managed_land, _ = init_crop(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    crop_cal.sdate .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = sum(Array(soil.NO3))
    nh4_before = sum(Array(soil.NH4))

    fertilizer!(crop_cal, managed_land, crop, soil, 1)
    crop.fphu .= 0.3f0
    fertilizer!(crop_cal, managed_land, crop, soil, 2)

    total_input = sum(Array(soil.NO3)) - no3_before +
                  sum(Array(soil.NH4)) - nh4_before
    @test total_input ≈ 10.0f0 atol = 1.0f-6
    @test Array(crop.nfertilizer)[1] == 0.0f0
end
