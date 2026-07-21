using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL-compatible potential and actual LAI" begin
    cells = 4096
    crop = init_crop(Float32, cells, CuArray)
    pet = init_pet(Float32, cells, CuArray)

    crop.state.phenology.is_growing .= Int32(1)
    crop.state.phenology.growing_days .= Int32(20)
    crop.state.phenology.senescence .= true
    crop.state.canopy.lai .= 0.1f0
    crop.state.canopy.lai_npp_deficit .= 0.3f0
    crop.state.carbon.biomass .= 10.0f0
    crop.state.carbon.leaf .= 1.0f0
    crop.state.carbon.root .= 2.0f0
    crop.state.carbon.pool .= 7.0f0
    crop.state.nitrogen.sufficiency .= 1.0f0
    crop.state.water.sufficiency .= 1.0f0

    carbon_allocation!(cft1, crop)
    synchronize()
    @test all(==(0.1f0), Array(crop.state.canopy.lai))
    @test all(==(0.3f0), Array(crop.state.canopy.lai_npp_deficit))
    @test all(iszero, Array(crop.auxiliary.canopy.actual_lai))

    pet.par .= 20.0f0
    apar_crop!(cft1, crop, pet)
    pet.eeq .= 2.0f0
    rain = CUDA.fill(5.0f0, cells)
    interception!(crop, cft1, pet.eeq, rain)
    synchronize()
    @test all(iszero, Array(crop.auxiliary.canopy.fpar))
    @test all(iszero, Array(crop.auxiliary.canopy.apar))
    @test all(iszero, Array(crop.auxiliary.canopy.canopy_wet))
    @test all(iszero, Array(crop.fluxes.water.interception))
end
