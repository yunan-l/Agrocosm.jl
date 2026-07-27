using Agrocosm
using Terrarium
using Test

# Discrete crop-management events: sowing establishes a seeded stand and resets the phenological clock;
# harvest exports the grain, returns the residue to the soil litter (mass-conserving over the root
# zone), and clears the stand. The residue return reuses the root-fraction distribution, so the crop
# carbon removed at harvest is fully accounted for as exported grain + soil residue.
@testset "Crop management (sowing / harvest)" begin
    Δz = 0.1
    grid = ColumnGrid(CPU(), UniformSpacing(Δz = Δz, N = 10))
    model = CropModel(grid, crop_pft("maize"); soil_hydrology = SoilHydrology(eltype(grid)))
    integrator = initialize(model; initializers = (temperature = 18.0,))
    state = integrator.state
    calendar = CropCalendar(Float64; sowing_day = 120, harvest_day = 280, residue_fraction = 0.25)
    scalar(field) = interior(field)[1, 1, 1]
    integral(field) = sum(interior(field)) * Δz

    @testset "management_time maps days to seconds" begin
        @test Agrocosm.management_time(120, Float64) == 120 * Terrarium.seconds_per_day(Float64)
    end

    @testset "sowing seeds the stand and resets the clock" begin
        set!(state.phenological_heat_units, 500.0)   # leftover from a previous season
        set!(state.crop_biomass, 0.0)
        sow!(integrator, calendar)
        @test scalar(state.crop_biomass) ≈ calendar.seed_carbon
        @test scalar(state.crop_nitrogen) ≈ calendar.seed_nitrogen
        @test scalar(state.phenological_heat_units) == 0.0
    end

    @testset "harvest exports grain, returns residue, clears the stand (mass-conserving)" begin
        # Establish a developed canopy, then a standing biomass to partition into organs.
        set!(state.phenological_heat_units, 0.6 * model.vegetation.phenology_dynamics.heat_unit_requirement)
        set!(state.crop_biomass, 0.8)
        set!(state.crop_nitrogen, 0.8 / 30)
        Terrarium.compute_auxiliary!(state, integrator.model)

        leaf_c = scalar(state.leaf_carbon)
        root_c = scalar(state.root_carbon)
        storage_c = scalar(state.storage_carbon)
        leaf_n = scalar(state.leaf_nitrogen)
        root_n = scalar(state.root_nitrogen)
        storage_n = scalar(state.storage_nitrogen)
        biomass_before = leaf_c + root_c + storage_c
        litter_before = integral(state.litter_carbon)
        litter_nitrogen_before = integral(state.litter_nitrogen)

        # Harvest index (grain fraction of aboveground) from the phenology + water status.
        fphu = clamp(scalar(state.phenology_heat_unit_fraction), 0.0, 1.0)
        hi = crop_harvest_index(calendar.harvest_index, fphu, 100.0 * scalar(state.soil_moisture_limiting_factor))
        aboveground_c = storage_c + leaf_c
        aboveground_n = storage_n + leaf_n

        yield = harvest!(integrator, calendar)

        # Stand cleared.
        @test scalar(state.crop_biomass) == 0.0
        @test scalar(state.crop_nitrogen) == 0.0
        @test scalar(state.phenological_heat_units) == 0.0

        # Grain export (HI × aboveground) and non-grain residue return ((1−HI)×aboveground + root).
        # The residue nitrogen enters the litter nitrogen pool (organic N), not the mineral pool.
        expected_yield = hi * aboveground_c
        residue_carbon = (1 - hi) * aboveground_c + root_c
        residue_nitrogen = (1 - hi) * aboveground_n + root_n
        @test yield ≈ expected_yield rtol = 1e-10
        @test integral(state.litter_carbon) - litter_before ≈ residue_carbon rtol = 1e-8
        @test integral(state.litter_nitrogen) - litter_nitrogen_before ≈ residue_nitrogen rtol = 1e-8

        # Carbon closes: everything removed from the crop is either exported (grain) or in the soil.
        @test yield + (integral(state.litter_carbon) - litter_before) ≈ biomass_before rtol = 1e-8
    end

    @testset "continuous fertilizer flux is windowed and split by the clock" begin
        # The initialized clock sits at t = 0 (day 0), so a window opening at day 0 is active and one
        # opening later is not — no clock mutation needed.
        rate = 1.0e-7
        active = CropFertilization(Float64; application_rate = rate, nitrate_fraction = 0.4,
            application_start_day = 0, application_end_day = 10)
        fertilize!(integrator, active)
        @test scalar(state.fertilizer_ammonium_flux) ≈ rate * 0.6
        @test scalar(state.fertilizer_nitrate_flux) ≈ rate * 0.4

        upcoming = CropFertilization(Float64; application_rate = rate, nitrate_fraction = 0.4,
            application_start_day = 5, application_end_day = 10)
        fertilize!(integrator, upcoming)
        @test scalar(state.fertilizer_ammonium_flux) == 0.0
        @test scalar(state.fertilizer_nitrate_flux) == 0.0
    end
end
