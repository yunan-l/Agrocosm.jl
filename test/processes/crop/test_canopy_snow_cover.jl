using Agrocosm
using Test

@testset "LPJmL crop snow cover suppresses absorbed PAR" begin
    for T in (Float32, Float64), pft in (cft1, cft3)
        crop = init_crop(T, 2, identity)
        pet = init_pet(T, 2, identity)
        crop.state.canopy.lai .= T(2)
        pet.par .= T(20)
        snow_height = T[0, 0.1]

        if pft === cft3
            apar_crop_maize!(pft, crop, pet, snow_height)
        else
            apar_crop!(pft, crop, pet, snow_height)
        end

        @test crop.auxiliary.canopy.fpar[1] > zero(T)
        @test crop.auxiliary.canopy.apar[1] > zero(T)
        @test crop.auxiliary.canopy.fpar[2] == zero(T)
        @test crop.auxiliary.canopy.apar[2] == zero(T)
    end
end
