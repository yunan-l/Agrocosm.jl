using Agrocosm
using Test

@testset "Fertilizer process has no implicit unlimited-N source" begin
    crop, crop_cal, managed_land, _ = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    no3_before = copy(soil.NO3)
    nh4_before = copy(soil.NH4)

    fertilizer!(crop_cal, managed_land, crop, soil, 1)

    @test soil.NO3 == no3_before
    @test soil.NH4 == nh4_before
end

@testset "Prescribed fertilizer is split and conserved" begin
    crop, crop_cal, managed_land, _ = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    crop_cal.sdate .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = sum(soil.NO3)
    nh4_before = sum(soil.NH4)

    fertilizer!(crop_cal, managed_land, crop, soil, 1)
    @test crop.nfertilizer[1] ≈ 8.0f0

    crop.fphu .= 0.3f0
    fertilizer!(crop_cal, managed_land, crop, soil, 2)

    no3_input = sum(soil.NO3) - no3_before
    nh4_input = sum(soil.NH4) - nh4_before
    @test no3_input ≈ 5.0f0 atol = 1.0f-6
    @test nh4_input ≈ 5.0f0 atol = 1.0f-6
    @test no3_input + nh4_input ≈ 10.0f0 atol = 1.0f-6
    @test crop.nfertilizer[1] == 0.0f0
end

@testset "Automatic-fertilizer mode disables prescribed inputs" begin
    crop, crop_cal, managed_land, _ = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    crop_cal.sdate .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = copy(soil.NO3)
    nh4_before = copy(soil.NH4)

    fertilizer!(crop_cal, managed_land, crop, soil, 1; enabled = false)

    @test soil.NO3 == no3_before
    @test soil.NH4 == nh4_before
    @test crop.nfertilizer[1] == 0.0f0
end
