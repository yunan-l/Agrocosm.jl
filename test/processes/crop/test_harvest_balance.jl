using Agrocosm
using Test

@testset "Harvest carbon and nitrogen routing is conservative" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    state = test_model_state(crop, soil; output)
    residue_fraction = Float32[0.6]
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]

    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 4.0f0
    crop.state.carbon.storage .= 5.0f0
    crop.state.carbon.pool .= 3.0f0
    crop.state.carbon.biomass .= 14.0f0
    crop.state.nitrogen.leaf .= 0.2f0
    crop.state.nitrogen.root .= 0.4f0
    crop.state.nitrogen.storage .= 0.5f0
    crop.state.nitrogen.pool .= 0.3f0
    crop.state.nitrogen.total .= 1.4f0
    crop.state.phenology.is_growing .= Int32(1)
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true

    carbon_before = crop.state.carbon.leaf[1] + crop.state.carbon.root[1] +
        crop.state.carbon.storage[1] + crop.state.carbon.pool[1]
    nitrogen_before = crop.state.nitrogen.total[1]
    harvest_crop!(state, state, output, residue_fraction, 150)

    carbon_export = crop.fluxes.carbon.harvest_export[1]
    nitrogen_export = crop.fluxes.nitrogen.harvest_export[1]
    carbon_residue = sum(soil.carbon.input)
    nitrogen_residue = sum(soil.nitrogen.input)
    @test crop.events.harvest[1] == 1
    @test crop.fluxes.carbon.yield[1] == 5.0f0
    @test carbon_export == 7.0f0
    @test carbon_residue == 7.0f0
    @test carbon_export + carbon_residue == carbon_before
    @test nitrogen_export == 0.7f0
    @test nitrogen_residue ≈ 0.7f0 atol = 1.0f-7
    @test nitrogen_export + nitrogen_residue ≈ nitrogen_before atol = 2.0f-7

    Agrocosm.route_harvest_residues!(state, state)
    @test sum(soil.carbon.litter) == carbon_residue
    @test sum(soil.nitrogen.litter) == nitrogen_residue
    @test crop.state.phenology.is_growing[1] == 0
end

@testset "Harvest-season diagnostics are accumulated and emitted annually" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    state = test_model_state(crop, soil; output)
    Agrocosm.prepare_output_block!(output, 1, 1)

    crop.events.sowing .= Int32(1)
    crop.state.phenology.is_growing .= Int32(1)
    crop.fluxes.carbon.gross_assimilation .= 4.0f0
    crop.auxiliary.canopy.actual_lai .= 2.0f0
    crop.auxiliary.stress.water_deficit .= 3.0f0
    crop.fluxes.water.interception .= 0.5f0
    crop.fluxes.water.transpiration_layer .= 0.1f0
    soil.water.evaporation .= 0.2f0
    soil.surface_litter.evaporation .= 0.3f0

    Agrocosm.accumulate_season_process_diagnostics!(output, state, state)
    crop.events.sowing .= Int32(0)
    Agrocosm.accumulate_season_process_diagnostics!(output, state, state)

    @test output.annual.active_gpp[1] == 8.0f0
    @test output.annual.active_lai_days[1] == 4.0f0
    @test output.annual.active_length[1] == 2.0f0
    @test output.annual.active_water_deficit[1] == 6.0f0
    @test output.annual.active_evapotranspiration[1] ≈ 4.6f0

    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.pool .= 3.0f0
    crop.state.carbon.storage .= 5.0f0
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true
    harvest_crop!(state, state, output, Float32[0.6], 365;
                  output_row = 1, annual_output_row = 1)

    @test output.crop.yield[1, 1] == 5.0f0
    @test output.crop.season_gpp[1, 1] == 8.0f0
    @test output.crop.season_lai_days[1, 1] == 4.0f0
    @test output.crop.season_length[1, 1] == 2.0f0
    @test output.crop.season_water_deficit[1, 1] == 6.0f0
    @test output.crop.season_evapotranspiration[1, 1] ≈ 4.6f0
    @test output.crop.harvest_aboveground_carbon[1, 1] == 10.0f0
end

@testset "Negative biomass terminates the failed crop conservatively" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    state = test_model_state(crop, soil; output)
    Agrocosm.prepare_output_block!(output, 1, 0)
    residue_fraction = Float32[0.6]
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]

    crop.state.phenology.is_growing .= Int32(1)
    crop.state.canopy.lai .= 1.0f0
    crop.auxiliary.canopy.actual_lai .= 1.0f0
    crop.state.carbon.biomass .= 1.0f0
    crop.state.carbon.leaf .= 0.2f0
    crop.state.carbon.root .= 0.3f0
    crop.state.carbon.pool .= 0.5f0
    crop.fluxes.carbon.respiration .= 2.0f0
    crop.state.nitrogen.total .= 0.4f0
    crop.state.nitrogen.leaf .= 0.1f0
    crop.state.nitrogen.root .= 0.1f0
    crop.state.nitrogen.pool .= 0.2f0

    carbon_allocation!(cft1, state)
    @test crop.fluxes.carbon.npp[1] == -2.0f0
    @test crop.events.harvest[1] == 1

    Agrocosm.terminate_failed_crop!(
        state, state, output, residue_fraction, 100; output_row = 1,
    )
    Agrocosm.route_harvest_residues!(state, state)

    @test crop.state.phenology.is_growing[1] == 0
    @test crop.events.harvest[1] == 1
    @test crop.state.carbon.biomass[1] == 0.0f0
    @test crop.state.carbon.leaf[1] == 0.0f0
    @test crop.state.carbon.root[1] == 0.0f0
    @test crop.state.carbon.pool[1] == 0.0f0
    @test crop.state.carbon.storage[1] == 0.0f0
    @test crop.state.nitrogen.total[1] == 0.0f0
    @test output.crop.growing_mask[1, 1] == 0
    @test output.calendar.harvest_event[1, 1] == 1
    @test crop.fluxes.carbon.harvest_export[1] ≈ -1.3f0
    @test all(>=(0), soil.carbon.input)
    @test all(>=(0), soil.carbon.litter)
    @test crop.fluxes.carbon.harvest_export[1] + sum(soil.carbon.input) ≈ -1.0f0
    @test crop.fluxes.nitrogen.harvest_export[1] + sum(soil.nitrogen.input) ≈ 0.4f0

    # An absent stand cannot continue producing a negative-NPP tail.
    carbon_allocation!(cft1, state)
    @test crop.fluxes.carbon.npp[1] == 0.0f0
end

@testset "Seed nitrogen is conserved when a crop fails on its sowing day" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    managed_land = init_managed_land(1, identity)
    state = test_model_state(crop, soil; managed_land, output)
    Agrocosm.prepare_output_block!(output, 1, 0)
    residue_fraction = Float32[0.6]
    crop.auxiliary.calendar.sowing_date .= Int32(100)
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]

    cultivate!(
        state, managed_land, state, 100;
        apply_prescribed_fertilizer = false,
    )
    crop.events.harvest .= Int32(1)
    Agrocosm.terminate_failed_crop!(
        state, state, output, residue_fraction, 100; output_row = 1,
    )
    Agrocosm.route_harvest_residues!(state, state)

    @test crop.state.phenology.is_growing[1] == 0
    @test crop.fluxes.nitrogen.harvest_export[1] + sum(soil.nitrogen.litter) ≈
          0.7f0 atol = 2.0f-7
end

function test_last_harvest_yield(device)
    @testset "Annual yield retains both sowing-year cohorts" begin
        for annual_row in (nothing, 1)
            crop = init_crop(4, device)
            soil = init_soil(4, soilparams.soildepth, device)
            output = init_output(4, device)
            state = test_model_state(crop, soil; output)
            Agrocosm.prepare_output_block!(output, 1, isnothing(annual_row) ? 0 : 1)
            residue_fraction = device(fill(0.6f0, 4))

            copyto!(crop.state.carbon.storage, Float32[11, 7, 0, 0])
            copyto!(crop.state.phenology.growing_days, Int32[120, 90, 0, 0])
            copyto!(crop.state.phenology.harvesting, Bool[true, true, false, false])
            crop.state.phenology.harvesting_previous .= false
            crop.state.phenology.is_growing .= Int32(1)
            harvest_crop!(state, state, output, residue_fraction, 90; output_row = 1)
            @test Array(output.annual.yield) == Float32[11, 7, 0, 0]
            @test Array(output.annual.yield_from_previous_sowing_year) == Float32[11, 0, 0, 0]
            @test Array(output.annual.yield_from_current_sowing_year) == Float32[0, 7, 0, 0]

            copyto!(crop.state.carbon.storage, Float32[5, 0, 9, 0])
            copyto!(crop.state.phenology.growing_days, Int32[200, 90, 100, 0])
            copyto!(crop.state.phenology.harvesting_previous, Bool[false, true, false, false])
            copyto!(crop.state.phenology.harvesting, Bool[true, true, true, false])
            crop.state.phenology.is_growing .= Int32(1)
            harvest_crop!(state, state, output, residue_fraction, 300; output_row = 1)
            @test Array(output.annual.yield) == Float32[5, 7, 9, 0]
            @test Array(output.annual.yield_from_previous_sowing_year) == Float32[11, 0, 0, 0]
            @test Array(output.annual.yield_from_current_sowing_year) == Float32[5, 7, 9, 0]
            @test Array(output.annual.harvest_date) == Int32[300, 90, 300, 0]

            crop.state.phenology.harvesting_previous .= true
            harvest_crop!(state, state, output, residue_fraction, 365;
                          output_row = 1, annual_output_row = annual_row)
            @test Array(output.crop.yield) == reshape(Float32[5, 7, 9, 0], 1, :)
            @test Array(output.crop.yield_from_previous_sowing_year) ==
                reshape(Float32[11, 0, 0, 0], 1, :)
            @test Array(output.crop.yield_from_current_sowing_year) ==
                reshape(Float32[5, 7, 9, 0], 1, :)
            @test Array(output.calendar.harvest_date) == reshape(Int32[300, 90, 300, 0], 1, :)
            @test all(iszero, Array(output.annual.yield))
            @test all(iszero, Array(output.annual.yield_from_previous_sowing_year))
            @test all(iszero, Array(output.annual.yield_from_current_sowing_year))
            @test all(iszero, Array(output.annual.harvest_date))
        end
    end

    @testset "A later failed harvest replaces the earlier yield" begin
        for failure_day in (300, 365)
            crop = init_crop(2, device)
            soil = init_soil(2, soilparams.soildepth, device)
            output = init_output(2, device)
            state = test_model_state(crop, soil; output)
            Agrocosm.prepare_output_block!(output, 1, 1)
            residue_fraction = device(fill(0.6f0, 2))

            copyto!(crop.state.carbon.storage, Float32[11, 13])
            crop.state.phenology.harvesting .= true
            crop.state.phenology.harvesting_previous .= false
            crop.state.phenology.is_growing .= Int32(1)
            harvest_crop!(state, state, output, residue_fraction, 90; output_row = 1)

            copyto!(crop.state.carbon.storage, Float32[3, -2])
            crop.state.phenology.harvesting .= false
            crop.state.phenology.is_growing .= Int32(1)
            if failure_day == 365
                harvest_crop!(state, state, output, residue_fraction, 365;
                              output_row = 1, annual_output_row = 1)
                @test Array(output.crop.yield) == reshape(Float32[11, 13], 1, :)
            end
            crop.events.harvest .= Int32(1)
            Agrocosm.terminate_failed_crop!(
                state, state, output, residue_fraction, failure_day;
                output_row = 1, annual_output_row = failure_day == 365 ? 1 : nothing,
            )
            if failure_day != 365
                @test Array(output.annual.yield) == Float32[3, -2]
                harvest_crop!(state, state, output, residue_fraction, 365;
                              output_row = 1, annual_output_row = 1)
            end
            @test Array(output.crop.yield) == reshape(Float32[3, 0], 1, :)
            @test Array(output.crop.yield_from_previous_sowing_year) ==
                zeros(Float32, 1, 2)
            @test Array(output.crop.yield_from_current_sowing_year) ==
                reshape(Float32[14, 13], 1, :)
            @test Array(output.calendar.harvest_date) == fill(Int32(failure_day), 1, 2)
            @test all(iszero, Array(output.annual.yield))
            @test all(iszero, Array(crop.state.carbon.storage))
        end
    end
end

test_last_harvest_yield(identity)
