using Agrocosm
using Test

@testset "Minimal crop checkpoint keeps only prognostic process state" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    managed_land = init_managed_land(1, identity)
    pet = init_pet(1, identity)

    restart = crop_restart_payload(crop)
    @test propertynames(crop.state) == (:phenology, :canopy, :carbon, :nitrogen, :water)
    @test :fphu ∉ propertynames(crop.state.phenology)
    @test :phenology_fraction ∉ propertynames(crop.state.canopy)
    @test :calendar ∉ propertynames(crop.state)
    @test restart.state === crop.state

    # The fertilizer split must derive fphu from prognostic husum and PHU,
    # rather than use a stale diagnostic cache.
    crop.state.nitrogen.pending_fertilizer .= 10.0f0
    crop.state.phenology.husum .= 300.0f0
    crop.auxiliary.phenology.phu .= 1000.0f0
    crop.auxiliary.phenology.fphu .= 0.0f0
    fertilizer!(crop, managed_land, soil, 2)
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] == 10.0f0
    @test crop.state.nitrogen.pending_fertilizer[1] == 0.0f0

    # Albedo receives the previous-day canopy state and reconstructs the
    # scalar phenology fraction locally as lai / laimax.
    crop.state.canopy.lai .= 0.4f0 * cft1.laimax
    crop.state.phenology.is_growing .= 1
    albedo!(cft1, crop, pet)
    expected_canopy = cft1.fpc * (0.4f0 * cft1.albedo_leaf + 0.6f0 * cft1.albedo_litter)
    @test crop.auxiliary.canopy.albedo[1] ≈ expected_canopy

    # Annual harvest records are output bookkeeping: harvesting changes no
    # prognostic calendar state, but retains report data until day 365.
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true
    crop.state.carbon.storage .= 5.0f0
    harvest_crop!(crop, soil, output, Float32[0.67], 100)
    @test output.annual.yield[1] == 5.0f0
    @test output.annual.harvest_date[1] == 100
end
