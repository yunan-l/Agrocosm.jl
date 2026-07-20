using Agrocosm
using Test

include("../../helpers/crop_lifecycle_fixture.jl")

@testset "LPJmL-style seasonal state is rebuilt at cultivation" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 2, identity)
        soil = init_soil(T, 2, soilparams.soildepth, identity)
        managed_land = init_managed_land(T, 2, identity)
        crop.calendar.sowing_date .= Int32[100, 101]

        # Emulate stale state from a completed previous crop season.
        for values in (
            crop.phenology.vdsum, crop.phenology.husum, crop.phenology.fphu,
            crop.canopy.lai_npp_deficit, crop.nitrogen.leaf,
            crop.nitrogen.root, crop.nitrogen.pool, crop.nitrogen.storage,
            crop.nitrogen.pending_fertilizer, crop.nitrogen.pending_manure,
            crop.nitrogen.stress_sum, crop.water.demand_sum,
            crop.water.supply_sum, crop.water.waterlogging_days,
        )
            values .= T(9)
        end
        crop.phenology.senescence .= true
        crop.phenology.senescence_previous .= true
        crop.phenology.growing_days .= Int32(99)
        crop.nitrogen.stress .= T(0.2)
        crop.water.stress .= T(0.3)
        crop.water.waterlogging_stress .= T(0.4)
        crop.carbon.yield .= T(7)
        soil_carbon_before = copy(soil.carbon.slow)

        cultivate!(
            crop, crop.calendar, managed_land, soil, 100;
            apply_prescribed_fertilizer = false,
            laimax = cft1.laimax,
        )

        @test crop.calendar.sowing_callback == Int32[1, 0]
        @test crop.phenology.is_growing[1] == 1
        @test crop.phenology.growing_days[1] == 0
        @test crop.phenology.vdsum[1] == zero(T)
        @test crop.phenology.husum[1] == zero(T)
        @test crop.phenology.fphu[1] == zero(T)
        @test !crop.phenology.senescence[1]
        @test !crop.phenology.senescence_previous[1]
        @test crop.canopy.flaimax[1] == T(0.000083)
        @test crop.canopy.lai[1] == T(0.000083) * T(cft1.laimax)
        @test crop.canopy.laimax_adjusted[1] == one(T)
        @test crop.nitrogen.total[1] == T(0.7)
        @test crop.nitrogen.leaf[1] == zero(T)
        @test crop.nitrogen.pending_fertilizer[1] == zero(T)
        @test crop.nitrogen.stress_sum[1] == zero(T)
        @test crop.nitrogen.stress[1] == one(T)
        @test crop.water.demand_sum[1] == zero(T)
        @test crop.water.supply_sum[1] == zero(T)
        @test crop.water.stress[1] == one(T)
        @test crop.water.waterlogging_days[1] == zero(T)
        @test crop.water.waterlogging_stress[1] == one(T)

        # Cultivation rebuilds only the new crop; annual output and soil memory persist.
        @test crop.carbon.yield == T[7, 7]
        @test soil.carbon.slow == soil_carbon_before
        @test crop.phenology.husum[2] == T(9)
        @test crop.nitrogen.pending_fertilizer[2] == T(9)
    end
end

@testset "Two-year prescribed crop lifecycle has one ordered event sequence per year" begin
    simulation = run_lifecycle_fixture()
    sowing_days = callback_days(simulation.output.calendar.sowing_callback)
    harvest_days = callback_days(simulation.output.calendar.harvest_callback)
    growing = vec(Array(simulation.output.crop.growing_mask))
    fphu = vec(Array(simulation.output.crop.fphu))

    @test sowing_days == [100, 465]
    @test harvest_days == [137, 502]
    @test length(sowing_days) == length(harvest_days) == 2
    for (sowing_day, harvest_day) in zip(sowing_days, harvest_days)
        @test sowing_day < harvest_day
        @test all(==(1), growing[sowing_day:(harvest_day - 1)])
        @test growing[harvest_day] == 0
        @test all(diff(fphu[sowing_day:(harvest_day - 1)]) .>= 0)
        @test fphu[sowing_day] > 0
        @test fphu[harvest_day - 1] == 1
    end
    @test all(==(0), growing[1:99])
    @test all(==(0), growing[137:464])
    @test all(==(0), growing[502:730])
    @test size(simulation.output.crop.yield) == (2, 1)

    # LPJmL deletes the crop PFT after harvest. Agrocosm retains allocated GPU
    # arrays, so every active-crop stock, flux, and seasonal accumulator must be
    # zero from the following day until the next sowing event.
    inactive_days = vcat(138:464, 503:730)
    for field in (
        :gpp, :npp, :lambda, :potential_vmax, :vmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu,
    )
        values = Array(getproperty(simulation.output.crop, field))
        @test all(iszero, values[inactive_days, :])
    end
    @test all(iszero, Array(simulation.output.crop.growing_mask)[inactive_days, :])
    for field in (:harvesting_mask, :sowing_callback, :harvest_callback)
        values = Array(getproperty(simulation.output.calendar, field))
        @test all(iszero, values[inactive_days, :])
    end

    # Final day 730 is fallow. Static templates and environmental diagnostics
    # are excluded: phu/winter_type, initial_organs, root_distribution, albedo,
    # root_zone_water, temperature responses, and neutral waterlogging stress.
    for (container, fields) in (
        (simulation.crop.phenology,
         (:vdsum, :husum, :fphu, :growing_days, :is_growing)),
        (simulation.crop.canopy,
         (:lai, :flaimax, :laimax_adjusted, :lai_npp_deficit,
          :phenology_fraction, :fpar, :apar)),
        (simulation.crop.carbon,
         (:biomass, :leaf, :root, :pool, :storage, :organs, :yield,
          :npp, :respiration)),
        (simulation.crop.nitrogen,
         (:total, :uptake, :auto_fertilizer, :leaf, :root, :pool, :storage,
          :demand_total, :demand_leaf, :pending_manure, :pending_fertilizer,
          :seed_input, :prescribed_manure_input,
          :prescribed_fertilizer_input, :harvest_export, :stress_sum,
          :stress, :deficit)),
        (simulation.crop.water,
         (:canopy_conductance, :transpiration, :canopy_wet, :interception,
          :transpiration_layer, :deficit, :demand_sum, :supply_sum, :stress,
          :waterlogging_days)),
        (simulation.crop.photosynthesis,
         (:gross_assimilation, :net_assimilation,
          :water_limited_assimilation, :leaf_respiration, :potential_vmax,
          :vmax, :nitrogen_limitation, :lambda)),
    )
        for field in fields
            @test all(iszero, Array(getproperty(container, field)))
        end
    end
    @test !simulation.crop.phenology.senescence[1]
    @test !simulation.crop.phenology.senescence_previous[1]
    @test simulation.crop.phenology.harvesting[1]
    @test simulation.crop.phenology.harvesting_previous[1]
    @test simulation.crop.water.waterlogging_stress[1] == 1.0f0
end

@testset "Cultivation, fertilizer, and tillage are single-trigger events" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    managed_land = init_managed_land(1, identity)
    crop.calendar.sowing_date .= Int32(100)
    managed_land.fertilizer .= 10.0f0
    soil.management.tillage_fraction .= Float32[
        0.05 0 0
        0.95 1 0
        0 0 1
    ]
    soil.carbon.litter[:, 1] .= Float32[100, 0, 0]
    carbon_before = sum(soil.carbon.litter)
    mineral_before = sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium)

    cultivate!(
        crop, crop.calendar, managed_land, soil, 99;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop.calendar)
    @test crop.calendar.sowing_callback[1] == 0
    @test sum(soil.carbon.litter) == carbon_before
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) == mineral_before

    cultivate!(
        crop, crop.calendar, managed_land, soil, 100;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop.calendar)
    carbon_after_sowing = copy(soil.carbon.litter)
    @test crop.calendar.sowing_callback[1] == 1
    @test carbon_after_sowing[1, 1] < 100.0f0
    @test sum(carbon_after_sowing) == carbon_before
    @test crop.nitrogen.prescribed_fertilizer_input[1] == 2.0f0
    @test crop.nitrogen.pending_fertilizer[1] == 8.0f0

    crop.phenology.fphu .= 0.3f0
    cultivate!(
        crop, crop.calendar, managed_land, soil, 101;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop.calendar)
    @test crop.calendar.sowing_callback[1] == 0
    @test soil.carbon.litter == carbon_after_sowing
    @test crop.nitrogen.prescribed_fertilizer_input[1] == 8.0f0
    @test crop.nitrogen.pending_fertilizer[1] == 0.0f0
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) ≈
        mineral_before + 10.0f0

    cultivate!(
        crop, crop.calendar, managed_land, soil, 102;
        apply_prescribed_fertilizer = true,
    )
    @test crop.nitrogen.prescribed_fertilizer_input[1] == 0.0f0
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) ≈
        mineral_before + 10.0f0
end
