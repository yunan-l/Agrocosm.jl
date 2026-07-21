using Agrocosm
using Test

include("../../helpers/crop_lifecycle_fixture.jl")

@testset "LPJmL-style seasonal state is rebuilt at cultivation" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 2, identity)
        soil = init_soil(T, 2, soilparams.soildepth, identity)
        managed_land = init_managed_land(T, 2, identity)
        crop.auxiliary.calendar.sowing_date .= Int32[100, 101]

        # Emulate stale state from a completed previous crop season.
        for values in (
            crop.state.phenology.vdsum, crop.state.phenology.husum,
            crop.state.canopy.lai_npp_deficit, crop.state.nitrogen.leaf,
            crop.state.nitrogen.root, crop.state.nitrogen.pool, crop.state.nitrogen.storage,
            crop.state.nitrogen.pending_fertilizer, crop.state.nitrogen.pending_manure,
            crop.state.nitrogen.stress_sum, crop.state.water.demand_sum,
            crop.state.water.supply_sum,
        )
            values .= T(9)
        end
        crop.state.phenology.senescence .= true
        crop.state.phenology.senescence_previous .= true
        crop.state.phenology.growing_days .= Int32(99)
        crop.state.nitrogen.sufficiency .= T(0.2)
        crop.state.water.sufficiency .= T(0.3)
        soil_carbon_before = copy(soil.carbon.slow)

        cultivate!(
        crop, managed_land, soil, 100;
            apply_prescribed_fertilizer = false,
            laimax = cft1.laimax,
        )

        @test crop.events.sowing == Int32[1, 0]
        @test crop.state.phenology.is_growing[1] == 1
        @test crop.state.phenology.growing_days[1] == 0
        @test crop.state.phenology.vdsum[1] == zero(T)
        @test crop.state.phenology.husum[1] == zero(T)
        @test crop.auxiliary.phenology.fphu[1] == zero(T)
        @test !crop.state.phenology.senescence[1]
        @test !crop.state.phenology.senescence_previous[1]
        @test crop.auxiliary.canopy.flaimax[1] == T(0.000083)
        @test crop.state.canopy.lai[1] == T(0.000083) * T(cft1.laimax)
        @test crop.state.canopy.laimax_adjusted[1] == one(T)
        @test crop.state.nitrogen.total[1] == T(0.7)
        @test crop.state.nitrogen.leaf[1] == zero(T)
        @test crop.state.nitrogen.pending_fertilizer[1] == zero(T)
        @test crop.state.nitrogen.stress_sum[1] == zero(T)
        @test crop.state.nitrogen.sufficiency[1] == one(T)
        @test crop.state.water.demand_sum[1] == zero(T)
        @test crop.state.water.supply_sum[1] == zero(T)
        @test crop.state.water.sufficiency[1] == one(T)

        # Cultivation rebuilds only the new crop; soil memory persists.
        @test soil.carbon.slow == soil_carbon_before
        @test crop.state.phenology.husum[2] == T(9)
        @test crop.state.nitrogen.pending_fertilizer[2] == T(9)
    end
end

@testset "Two-year prescribed crop lifecycle has one ordered event sequence per year" begin
    simulation = run_lifecycle_fixture()
    sowing_days = event_days(simulation.output.calendar.sowing_event)
    harvest_days = event_days(simulation.output.calendar.harvest_event)
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
    @test all(>(0), Array(simulation.output.crop.yield))

    # LPJmL deletes the crop PFT after harvest. Agrocosm retains allocated GPU
    # arrays, so every active-crop stock, flux, and seasonal accumulator must be
    # zero from the following day until the next sowing event.
    inactive_days = vcat(138:464, 503:730)
    for field in (
        :gpp, :npp, :lambda, :potential_vcmax, :vcmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu,
    )
        values = Array(getproperty(simulation.output.crop, field))
        @test all(iszero, values[inactive_days, :])
    end
    @test all(iszero, Array(simulation.output.crop.growing_mask)[inactive_days, :])
    for field in (:harvesting_mask, :sowing_event, :harvest_event)
        values = Array(getproperty(simulation.output.calendar, field))
        @test all(iszero, values[inactive_days, :])
    end

    # Final day 730 is fallow. Static templates and environmental diagnostics
    # are excluded: phu/winter_type, root_distribution, albedo and temperature responses.
    for (container, fields) in (
        (simulation.crop.state.phenology,
         (:vdsum, :husum, :growing_days, :is_growing)),
        (simulation.crop.state.canopy,
         (:lai, :laimax_adjusted, :lai_npp_deficit)),
        (simulation.crop.auxiliary.phenology, (:fphu,)),
        (simulation.crop.auxiliary.canopy,
         (:actual_lai, :flaimax, :fpar, :apar, :canopy_conductance, :canopy_wet)),
        (simulation.crop.state.carbon,
         (:biomass, :leaf, :root, :pool, :storage)),
        (simulation.crop.fluxes.carbon,
         (:yield, :harvest_export, :npp, :respiration, :gross_assimilation, :net_assimilation,
          :water_limited_assimilation, :leaf_respiration)),
        (simulation.crop.state.nitrogen,
         (:total, :leaf, :root, :pool, :storage, :pending_manure,
          :pending_fertilizer, :stress_sum)),
        (simulation.crop.fluxes.nitrogen,
         (:uptake, :auto_fertilizer, :seed_input, :prescribed_manure_input,
          :prescribed_fertilizer_input, :harvest_export)),
        (simulation.crop.state.water,
         (:demand_sum, :supply_sum)),
        (simulation.crop.fluxes.water,
         (:interception, :transpiration_layer)),
        (simulation.crop.auxiliary.stress,
         (:nitrogen_demand_total, :nitrogen_demand_leaf,
          :nitrogen_deficit, :water_deficit)),
        (simulation.crop.auxiliary.photosynthesis,
         (:potential_vcmax, :vcmax, :nitrogen_limitation, :lambda)),
    )
        for field in fields
            @test all(iszero, Array(getproperty(container, field)))
        end
    end
    @test !simulation.crop.state.phenology.senescence[1]
    @test !simulation.crop.state.phenology.senescence_previous[1]
    @test !simulation.crop.state.phenology.harvesting[1]
    @test !simulation.crop.state.phenology.harvesting_previous[1]
    @test simulation.output.annual.yield[1] == 0.0f0
    @test simulation.crop.state.nitrogen.sufficiency[1] == 1.0f0
    @test simulation.crop.state.water.sufficiency[1] == 1.0f0
end

@testset "Cultivation, fertilizer, and tillage are single-trigger events" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    managed_land = init_managed_land(1, identity)
    crop.auxiliary.calendar.sowing_date .= Int32(100)
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
        crop, managed_land, soil, 99;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop)
    @test crop.events.sowing[1] == 0
    @test sum(soil.carbon.litter) == carbon_before
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) == mineral_before

    cultivate!(
        crop, managed_land, soil, 100;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop)
    carbon_after_sowing = copy(soil.carbon.litter)
    @test crop.events.sowing[1] == 1
    @test carbon_after_sowing[1, 1] < 100.0f0
    @test sum(carbon_after_sowing) == carbon_before
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] == 2.0f0
    @test crop.state.nitrogen.pending_fertilizer[1] == 8.0f0

    crop.auxiliary.phenology.phu .= 1.0f0
    crop.state.phenology.husum .= 0.3f0
    cultivate!(
        crop, managed_land, soil, 101;
        apply_prescribed_fertilizer = true,
    )
    litter_tillage!(soil, crop)
    @test crop.events.sowing[1] == 0
    @test soil.carbon.litter == carbon_after_sowing
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] == 8.0f0
    @test crop.state.nitrogen.pending_fertilizer[1] == 0.0f0
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) ≈
        mineral_before + 10.0f0

    cultivate!(
        crop, managed_land, soil, 102;
        apply_prescribed_fertilizer = true,
    )
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] == 0.0f0
    @test sum(soil.nitrogen.nitrate) + sum(soil.nitrogen.ammonium) ≈
        mineral_before + 10.0f0
end
