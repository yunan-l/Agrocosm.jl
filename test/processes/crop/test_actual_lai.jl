using Agrocosm
using Test

@testset "LPJmL-compatible potential and actual LAI" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 1, identity)
        pet = init_pet(T, 1, identity)

        crop.state.phenology.is_growing .= Int32(1)
        crop.state.phenology.growing_days .= Int32(20)
        crop.state.phenology.senescence .= true
        crop.state.canopy.lai .= T(0.1)
        crop.state.canopy.lai_npp_deficit .= T(0.3)
        crop.state.carbon.biomass .= T(10)
        crop.state.carbon.leaf .= T(1)
        crop.state.carbon.root .= T(2)
        crop.state.carbon.pool .= T(7)
        crop.state.nitrogen.sufficiency .= one(T)
        crop.state.water.sufficiency .= one(T)

        carbon_allocation!(cft1, crop)

        # LPJmL retains potential phenological LAI and clips only actual LAI.
        @test only(crop.state.canopy.lai) == T(0.1)
        @test only(crop.state.canopy.lai_npp_deficit) == T(0.3)
        @test only(crop.auxiliary.canopy.actual_lai) == zero(T)

        pet.par .= T(20)
        apar_crop!(cft1, crop, pet)
        @test only(crop.auxiliary.canopy.fpar) == zero(T)
        @test only(crop.auxiliary.canopy.apar) == zero(T)

        pet.eeq .= T(2)
        rain = fill(T(5), 1)
        interception!(crop, cft1, pet.eeq, rain)
        @test only(crop.auxiliary.canopy.canopy_wet) == zero(T)
        @test only(crop.fluxes.water.interception) == zero(T)
    end
end
