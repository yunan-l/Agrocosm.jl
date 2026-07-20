using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA fused daily climate and crop kernels" begin
    cells = 4096
    days = 3
    temperature = reshape(Float32.(range(-5, 30; length = days * cells)), days, cells)
    precipitation = reshape(Float32.(range(0, 12; length = days * cells)), days, cells)
    shortwave = reshape(Float32.(range(0, 350; length = days * cells)), days, cells)
    longwave = reshape(Float32.(range(-100, 20; length = days * cells)), days, cells)
    wind = reshape(Float32.(range(1, 6; length = days * cells)), days, cells)
    climate_cpu = (
        temp = temperature, prec = precipitation, sw = shortwave,
        lw = longwave, wind = wind, co2 = Float32[400],
    )
    climate_gpu = (
        temp = CuArray(temperature), prec = CuArray(precipitation),
        sw = CuArray(shortwave), lw = CuArray(longwave),
        wind = CuArray(wind), co2 = CuArray(Float32[400]),
    )
    weather_reference = init_weather(cells, identity)
    weather_gpu = init_weather(cells, CuArray)
    Agrocosm.readclimate_reference!(climate_cpu, weather_reference, 2)
    readclimate!(climate_gpu, weather_gpu, 2)
    synchronize()
    @test Array(weather_gpu.temp) == weather_reference.temp
    @test Array(weather_gpu.wind) == weather_reference.wind

    daily_co2 = reshape(Float32.(range(390, 430; length = days * cells)), days, cells)
    daily_climate_cpu = merge(climate_cpu, (co2 = daily_co2,))
    daily_climate_gpu = merge(climate_gpu, (co2 = CuArray(daily_co2),))
    reference_co2 = Agrocosm.readclimate_reference!(
        daily_climate_cpu, weather_reference, 2,
    )
    gpu_co2 = readclimate!(daily_climate_gpu, weather_gpu, 2)
    synchronize()
    @test reference_co2 === weather_reference.daily_co2
    @test gpu_co2 === weather_gpu.daily_co2
    @test Array(gpu_co2) == reference_co2

    climate_bytes = CUDA.@allocated begin
        readclimate!(climate_gpu, weather_gpu, 2)
        synchronize()
    end

    crop_reference = init_crop(cells, identity)
    crop_gpu = init_crop(cells, CuArray)
    pet_reference = init_pet(cells, identity)
    pet_gpu = init_pet(cells, CuArray)
    lai = Float32.(range(0, 7; length = cells))
    phenology = Float32.(range(0, 1; length = cells))
    growing = Int32.(mod.(1:cells, 2))
    par = Float32.(range(0, 25; length = cells))
    crop_reference.state.canopy.lai .= lai
    crop_reference.auxiliary.canopy.phenology_fraction .= phenology
    crop_reference.state.phenology.is_growing .= growing
    pet_reference.par .= par
    crop_gpu.state.canopy.lai .= CuArray(lai)
    crop_gpu.auxiliary.canopy.phenology_fraction .= CuArray(phenology)
    crop_gpu.state.phenology.is_growing .= CuArray(growing)
    pet_gpu.par .= CuArray(par)
    Agrocosm.albedo_reference!(cft1, crop_reference, pet_reference)
    Agrocosm.apar_crop_reference!(cft1, crop_reference, pet_reference)
    albedo!(cft1, crop_gpu, pet_gpu)
    apar_crop!(cft1, crop_gpu, pet_gpu)
    synchronize()
    @test Array(pet_gpu.albedo) ≈ pet_reference.albedo rtol = 3.0f-6
    @test Array(crop_gpu.auxiliary.canopy.apar) ≈ crop_reference.auxiliary.canopy.apar rtol = 3.0f-6
    canopy_bytes = CUDA.@allocated begin
        albedo!(cft1, crop_gpu, pet_gpu)
        apar_crop!(cft1, crop_gpu, pet_gpu)
        synchronize()
    end

    land_reference = init_managed_land(cells, identity)
    land_gpu = init_managed_land(cells, CuArray)
    soil_reference = init_soil(cells, soilparams.soildepth, identity)
    soil_gpu = init_soil(cells, soilparams.soildepth, CuArray)
    sowing_dates = fill(Int32(101), cells)
    sowing_dates[1:2:end] .= Int32(100)
    crop_reference.state.calendar.sowing_date .= sowing_dates
    crop_gpu.state.calendar.sowing_date .= CuArray(sowing_dates)
    Agrocosm.cultivate_reference!(
        crop_reference, land_reference, soil_reference, 100;
        apply_prescribed_fertilizer = false,
    )
    cultivate!(
        crop_gpu, land_gpu, soil_gpu, 100;
        apply_prescribed_fertilizer = false,
    )
    synchronize()
    @test Array(crop_gpu.events.sowing) == crop_reference.events.sowing
    @test Array(crop_gpu.state.carbon.biomass) == crop_reference.state.carbon.biomass
    cultivation_bytes = CUDA.@allocated begin
        cultivate!(
        crop_gpu, land_gpu, soil_gpu, 100;
            apply_prescribed_fertilizer = false,
        )
        synchronize()
    end

    pet_gpu.daylength .= 12.0f0
    pet_gpu.eeq .= 4.0f0
    crop_gpu.auxiliary.canopy.fpar .= 0.6f0
    crop_gpu.auxiliary.canopy.canopy_wet .= 0.1f0
    crop_gpu.auxiliary.stress.root_distribution .= CuArray(Float32[0.35, 0.25, 0.18, 0.13, 0.09])
    soil_gpu.water.relative_content .= 0.55f0
    soil_gpu.water.holding_capacity_storage .= 100.0f0
    assimilation_gpu = CUDA.fill(8.0f0, cells)
    co2_gpu = CuArray(Float32[40])
    transpiration!(assimilation_gpu, cft1, crop_gpu, pet_gpu, soil_gpu, co2_gpu)
    synchronize()
    transpiration_bytes = CUDA.@allocated begin
        transpiration!(assimilation_gpu, cft1, crop_gpu, pet_gpu, soil_gpu, co2_gpu)
        synchronize()
    end

    @test climate_bytes == 0
    @test canopy_bytes == 0
    @test cultivation_bytes == 0
    @test transpiration_bytes == 0
end
