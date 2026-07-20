using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA state initialization" begin
    cell_size = 2
    # CUDA.jl does not guarantee that the UnionAll `CuArray` constructor is
    # fully inferred across Julia/CUDA.jl versions. Validate the concrete GPU
    # state fields below instead of testing CUDA.jl's internal inference.
    crop = init_crop(cell_size, CuArray)
    managed_land = init_managed_land(cell_size, CuArray)
    soil = init_soil(cell_size, soilparams.soildepth, CuArray)
    weather = init_weather(cell_size, CuArray)
    output = init_output(cell_size, CuArray)
    nitrogen_balance = init_nitrogen_balance(3, cell_size, CuArray)

    @test crop.phenology.phu isa CuArray{Float32, 1}
    @test crop.phenology.is_growing isa CuArray{Int32, 1}
    @test all(Array(crop.phenology.is_growing) .== Int32(0))
    @test crop.canopy.lai isa CuArray{Float32, 1}
    @test crop.carbon.organs isa CuArray{Float32, 2}
    @test crop.nitrogen.total isa CuArray{Float32, 1}
    @test crop.nitrogen.seed_input isa CuArray{Float32, 1}
    @test crop.nitrogen.harvest_export isa CuArray{Float32, 1}
    @test crop.water.transpiration_layer isa CuArray{Float32, 2}
    @test crop.calendar.sowing_date isa CuArray{Int32, 1}
    @test managed_land.latitude isa CuArray{Float32, 1}
    @test crop.photosynthesis.gross_assimilation isa CuArray{Float32, 1}
    @test soil.properties.sand_fraction isa CuArray{Float32, 2}
    @test soil.water.storage isa CuArray{Float32, 2}
    @test soil.water.ice_storage isa CuArray{Float32, 2}
    @test soil.water.wilting_ice_fraction isa CuArray{Float32, 2}
    @test soil.water.available_ice_storage isa CuArray{Float32, 2}
    @test soil.water.free_ice_storage isa CuArray{Float32, 2}
    @test soil.thermal.temperature isa CuArray{Float32, 2}
    @test soil.thermal.enthalpy isa CuArray{Float32, 2}
    @test soil.thermal.frozen_fraction isa CuArray{Float32, 2}
    @test soil.thermal.initialized isa CuArray{Bool, 1}
    @test soil.carbon.litter isa CuArray{Float32, 2}
    @test soil.carbon.litter_to_fast isa CuArray{Float32, 2}
    @test soil.carbon.litter_to_slow isa CuArray{Float32, 2}
    @test soil.nitrogen.nitrate isa CuArray{Float32, 2}
    @test soil.nitrogen.litter_to_fast isa CuArray{Float32, 2}
    @test soil.nitrogen.litter_to_slow isa CuArray{Float32, 2}
    @test soil.nitrogen.leaching isa CuArray{Float32, 1}
    @test soil.decomposition.response isa CuArray{Float32, 2}
    @test soil.management.tillage_fraction isa CuArray{Float32, 2}
    @test soil.surface_litter.water_storage isa CuArray{Float32, 1}
    @test soil.snow.pack isa CuArray{Float32, 1}
    @test weather.temp isa CuArray{Float32, 1}
    @test output.crop.npp isa CuArray{Float32, 2}
    @test output.crop.potential_vmax isa CuArray{Float32, 2}
    @test output.crop.nitrogen_limitation isa CuArray{Float32, 2}
    @test output.soil.water_storage isa CuArray{Float32, 2}
    @test output.climate.temperature isa CuArray{Float32, 2}
    @test output.calendar.harvest_date isa CuArray{Int32, 2}
    @test nitrogen_balance.residual isa CuArray{Float32, 2}

    @test size(crop.carbon.organs) == (4, cell_size)
    @test size(crop.water.transpiration_layer) == (5, cell_size)
    @test size(nitrogen_balance.residual) == (3, cell_size)
    @test all(Array(crop.canopy.lai) .== 0.0f0)

    u0 = (
        soil_NO3 = CUDA.fill(9000.0f0, 5, cell_size),
        soil_NH4 = CUDA.fill(8000.0f0, 5, cell_size),
    )
    soil.nitrogen.slow .= 100.0f0
    Agrocosm.initialize_soil_mineral_nitrogen!(soil, u0, :lpjml_initsoil)
    @test all(Array(soil.nitrogen.nitrate) .== 1.0f0)
    @test all(Array(soil.nitrogen.ammonium) .== 1.0f0)

    Agrocosm.initialize_soil_c_shift!(soil, (u0 = nothing,), :lpjml_initsoil)
    expected_shift = Float32[0.55, 0.1125, 0.1125, 0.1125, 0.1125]
    @test Array(soil.carbon.shift_fast) == repeat(expected_shift, 1, cell_size)
    @test Array(soil.carbon.shift_slow) == repeat(expected_shift, 1, cell_size)
    @test Array(soil.nitrogen.shift_fast) == Array(soil.carbon.shift_fast)
    @test Array(soil.nitrogen.shift_slow) == Array(soil.carbon.shift_slow)
end
