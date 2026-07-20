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
    @test length(soil.surface_litter.water_storage) == cell_size
    @test length(weather.temp) == cell_size
    @test output.crop isa CropOutput
    @test output.soil isa SoilOutput
    @test output.climate isa ClimateOutput
    @test output.calendar isa CalendarOutput
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
