using Agrocosm
using Test

@testset "Prescribed PHU maturity follows LPJmL daily event order" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    managed_land = init_managed_land(1, identity)
    crop.calendar.sowing_date .= Int32(100)
    crop.phenology.phu .= 543.0f0
    temperature = Float32[15]
    daylength = Float32[12]
    vernalization_requirement = Float32[0]
    senescence_days = Int[]
    harvest_days = Int[]

    for day in 100:137
        cultivate!(
            crop, crop.calendar, managed_land, soil, day;
            apply_prescribed_fertilizer = false,
            laimax = cft1.laimax,
        )
        phenology_crop!(
            crop, vernalization_requirement, cft1, temperature, daylength,
        )
        harvest_crop!(crop.calendar, crop, soil, output, Float32[0.67], day)
        crop.phenology.senescence[1] && push!(senescence_days, day)
        crop.calendar.harvest_callback[1] == 1 && push!(harvest_days, day)
    end

    senescence_increment = ceil(Int, cft1.fphusen * 543 / 15)
    # Senescence starts on the first daily increment crossing fphusen. Thirty-
    # seven increments reach PHU; LPJmL harvests from the prior husum on day 38.
    @test first(senescence_days) == 100 + senescence_increment - 1
    @test crop.phenology.fphu[1] == 1.0f0
    @test harvest_days == [137]
    @test crop.calendar.harvest_date[1] == 137
    @test crop.phenology.growing_days[1] == 38
    @test crop.phenology.is_growing[1] == 0
end
