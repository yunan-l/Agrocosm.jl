using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA prescribed fertilizer is conserved" begin
    crop = init_crop(1, CuArray)
    managed_land = init_managed_land(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = sum(Array(soil.nitrogen.nitrate))
    nh4_before = sum(Array(soil.nitrogen.ammonium))

    fertilizer!(state, managed_land, state, 1)
    @test Array(crop.fluxes.nitrogen.prescribed_fertilizer_input)[1] ≈ 2.0f0
    crop.auxiliary.phenology.phu .= 1.0f0
    crop.state.phenology.husum .= 0.3f0
    fertilizer!(state, managed_land, state, 2)
    @test Array(crop.fluxes.nitrogen.prescribed_fertilizer_input)[1] ≈ 8.0f0

    total_input = sum(Array(soil.nitrogen.nitrate)) - no3_before +
                  sum(Array(soil.nitrogen.ammonium)) - nh4_before
    @test total_input ≈ 10.0f0 atol = 1.0f-6
    @test Array(crop.state.nitrogen.pending_fertilizer)[1] == 0.0f0
end

@testset "CUDA LPJmL second manure dose remains on the surface" begin
    crop = init_crop(1, CuArray)
    managed_land = init_managed_land(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.manure .= 10.0f0

    fertilizer!(state, managed_land, state, 1; manure = true)
    incorporated_after_sowing = Array(soil.nitrogen.litter[2, :])
    crop.auxiliary.phenology.phu .= 1.0f0
    crop.state.phenology.husum .= 0.3f0
    fertilizer!(
        state, managed_land, state, 2;
        manure = true, surface_second_manure = true,
    )

    @test Array(soil.nitrogen.litter[1, :])[1] > 0.0f0
    @test Array(soil.carbon.litter[1, :])[1] > 0.0f0
    @test Array(soil.nitrogen.litter[2, :]) == incorporated_after_sowing
    @test Array(crop.state.nitrogen.pending_manure)[1] == 0.0f0
end
