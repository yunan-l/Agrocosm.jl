using Agrocosm
using Test

@testset "CPU state initialization" begin
    cell_size = 2
    crop = @inferred init_crop(cell_size, identity)
    managed_land = @inferred init_managed_land(cell_size, identity)
    soil = @inferred init_soil(cell_size, soilparams.soildepth, identity)
    weather = @inferred init_weather(cell_size, identity)
    pet = @inferred init_pet(cell_size, identity)
    climbuf = @inferred init_climbuf(cell_size, identity)
    output = @inferred init_output(cell_size, identity)

    @test size(crop.carbon.organs) == (4, cell_size)
    @test size(crop.water.transpiration_layer) == (5, cell_size)
    @test length(crop.calendar.sowing_date) == cell_size
    @test length(managed_land.latitude) == cell_size
    @test length(crop.photosynthesis.gross_assimilation) == cell_size
    @test length(pet.daylength) == cell_size
    @test size(climbuf.temp) == (31, cell_size)
    @test eltype(crop.canopy.lai) == Float32
    @test all(iszero, crop.canopy.lai)
    @test all(iszero, crop.photosynthesis.gross_assimilation)
    @test crop.phenology isa CropPhenology
    @test crop.canopy isa CropCanopy
    @test crop.carbon isa CropCarbon
    @test crop.nitrogen isa CropNitrogen
    @test crop.water isa CropWater
    @test crop.calendar isa CropCalendar
    @test crop.photosynthesis isa CropPhotosynthesis
    @test soil.properties isa SoilProperties
    @test soil.water isa SoilWater
    @test soil.thermal isa SoilThermal
    @test length(soil.thermal.initialized) == cell_size
    @test all(.!soil.thermal.initialized)
    @test size(soil.thermal.enthalpy) == (5, cell_size)
    @test size(soil.thermal.frozen_fraction) == (5, cell_size)
    @test soil.carbon isa SoilCarbon
    @test soil.nitrogen isa SoilNitrogen
    @test soil.decomposition isa SoilDecomposition
    @test size(soil.decomposition.litter_response) == (3, cell_size)
    @test soil.management isa SoilManagement
    @test soil.surface_litter isa SoilSurfaceLitter
    @test soil.snow isa SoilSnow
    @test size(soil.water.storage) == (5, cell_size)
    @test size(soil.water.ice_storage) == (5, cell_size)
    @test size(soil.water.wilting_ice_fraction) == (5, cell_size)
    @test size(soil.water.available_ice_storage) == (5, cell_size)
    @test size(soil.water.free_ice_storage) == (5, cell_size)
    @test size(soil.carbon.litter) == (3, cell_size)
    @test size(soil.carbon.litter_to_fast) == (5, cell_size)
    @test size(soil.carbon.litter_to_slow) == (5, cell_size)
    @test size(soil.nitrogen.litter_to_fast) == (5, cell_size)
    @test size(soil.nitrogen.litter_to_slow) == (5, cell_size)
    @test length(soil.surface_litter.water_storage) == cell_size
    @test length(weather.temp) == cell_size
    @test output.crop isa CropOutput
    @test output.soil isa SoilOutput
    @test output.climate isa ClimateOutput
    @test output.calendar isa CalendarOutput
end

@testset "LPJmL c_shift initialization strategies" begin
    cells = 2
    soil = init_soil(cells, soilparams.soildepth, identity)
    model_state = (u0 = nothing,)

    Agrocosm.initialize_soil_c_shift!(soil, model_state, :lpjml_initsoil)
    expected = Float32[0.55, 0.1125, 0.1125, 0.1125, 0.1125]
    @test soil.carbon.shift_fast == repeat(expected, 1, cells)
    @test soil.carbon.shift_slow == repeat(expected, 1, cells)
    @test soil.nitrogen.shift_fast == soil.carbon.shift_fast
    @test soil.nitrogen.shift_slow == soil.carbon.shift_slow
    @test vec(sum(soil.carbon.shift_fast; dims = 1)) ≈ ones(Float32, cells)

    restart_fast = repeat(Float32[0.4, 0.25, 0.15, 0.1, 0.1], 1, cells)
    restart_slow = repeat(Float32[0.5, 0.2, 0.15, 0.1, 0.05], 1, cells)
    restart_state = (c_shift_fast = restart_fast, c_shift_slow = restart_slow)
    Agrocosm.initialize_soil_c_shift!(soil, restart_state, :restart)
    @test soil.carbon.shift_fast == restart_fast
    @test soil.carbon.shift_slow == restart_slow
    @test soil.nitrogen.shift_fast == restart_fast
    @test soil.nitrogen.shift_slow == restart_slow

    @test_throws ArgumentError Agrocosm.initialize_soil_c_shift!(
        soil, model_state, :restart,
    )
    @test_throws ArgumentError Agrocosm.initialize_soil_c_shift!(
        soil, model_state, :unknown,
    )
end

@testset "Initial data loader makes c_shift optional" begin
    cells = 2
    layers = 5
    u0 = (
        swc = fill(100.0f0, layers, cells),
        litc = fill(1.0f0, 3, cells),
        fastc = fill(10.0f0, layers, cells),
        slowc = fill(100.0f0, layers, cells),
        litn = fill(0.1f0, 3, cells),
        fastn = fill(1.0f0, layers, cells),
        slown = fill(10.0f0, layers, cells),
    )
    initial_without_shift = (u0 = u0,)
    data_without_shift = (
        latitude = Float32[45, 46],
        crop = (
            sdate = Int32[100, 101],
            phu = Float32[1200, 1250],
            manure = zeros(Float32, cells),
            fertilizer = fill(50.0f0, cells),
            residuefrac = fill(0.3f0, cells),
        ),
        soilparam = (
            soilph = fill(6.5f0, cells),
            w_sat = fill(0.45f0, layers, cells),
            sand = fill(0.4f0, cells),
            clay = fill(0.2f0, cells),
            tdiff_0 = fill(0.2f0, cells),
            tdiff_15 = fill(0.5f0, cells),
            soildepth = Float32[200, 300, 500, 700, 1300],
        ),
        initialLPJmL = initial_without_shift,
    )

    loaded = InitialDataLoader(data_without_shift, [1, 2], identity)
    @test !hasproperty(loaded.ModelState, :c_shift_fast)
    @test !hasproperty(loaded.ModelState, :c_shift_slow)

    fast_shift = repeat(Float32[0.4, 0.25, 0.15, 0.1, 0.1], 1, cells)
    slow_shift = repeat(Float32[0.5, 0.2, 0.15, 0.1, 0.05], 1, cells)
    data_with_shift = merge(data_without_shift, (
        initialLPJmL = merge(initial_without_shift, (
            c_shift_fast = fast_shift,
            c_shift_slow = slow_shift,
        )),
    ))
    restart = InitialDataLoader(
        data_with_shift, [1, 2], identity; load_c_shift_restart = true,
    )
    @test restart.ModelState.c_shift_fast == fast_shift
    @test restart.ModelState.c_shift_slow == slow_shift
end

@testset "LPJmL mineral-N initialization strategies" begin
    slow_n = reshape(Float32[100, 200, 300, 400, 500], 5, 1)
    u0 = (
        soil_NO3 = fill(9000.0f0, 5, 1),
        soil_NH4 = fill(8000.0f0, 5, 1),
    )
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.nitrogen.slow .= slow_n

    Agrocosm.initialize_soil_mineral_nitrogen!(soil, u0, :restart)
    @test soil.nitrogen.nitrate == u0.soil_NO3
    @test soil.nitrogen.ammonium == u0.soil_NH4

    Agrocosm.initialize_soil_mineral_nitrogen!(soil, u0, :lpjml_initsoil)
    @test soil.nitrogen.nitrate ≈ slow_n ./ 100
    @test soil.nitrogen.ammonium ≈ slow_n ./ 100
    @test_throws ArgumentError Agrocosm.initialize_soil_mineral_nitrogen!(
        soil, u0, :unknown,
    )
    @test_throws ArgumentError Agrocosm.initialize_soil_mineral_nitrogen!(
        soil, (slown = slow_n,), :restart,
    )
end
