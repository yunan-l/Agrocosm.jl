using Agrocosm
using Test

@testset "Agrocosm.jl" begin
    @testset "CPU state initialization" begin
        cell_size = 2
        crop, crop_cal, managed_land, photos = init_crop(cell_size, identity)
        pet = init_pet(cell_size, identity)
        climbuf = init_climbuf(cell_size, identity)

        @test size(crop.vegc) == (4, cell_size)
        @test size(crop.trans_layer) == (5, cell_size)
        @test length(crop_cal.sdate) == cell_size
        @test length(managed_land.latitude) == cell_size
        @test length(photos.agd) == cell_size
        @test length(pet.daylength) == cell_size
        @test size(climbuf.temp) == (31, cell_size)
        @test eltype(crop.lai) == Float32
        @test all(iszero, crop.lai)
        @test all(iszero, photos.agd)
    end

    @testset "C3 photosynthesis CPU smoke test" begin
        _, _, _, photos = init_crop(1, identity)
        photos.tstress .= 1.0f0

        photosynthesis_C3!(
            cft1,
            photos,
            Float32[10.0],
            Float32[12.0],
            Float32[20.0],
            Float32[40.0];
            comp_vmax = true,
        )

        @test photos.lambda == Float32[0.8]
        @test all(isfinite, photos.vmax)
        @test all(isfinite, photos.agd)
        @test all(isfinite, photos.rd)
        @test all(isfinite, photos.adt)
        @test all(isfinite, photos.adtmm)
        @test all(photos.vmax .>= 0.0f0)
        @test all(photos.agd .>= 0.0f0)
        @test all(photos.adt .>= 0.0f0)

        photos.tstress .= 0.0f0
        photosynthesis_C3!(
            cft1,
            photos,
            Float32[10.0],
            Float32[12.0],
            Float32[20.0],
            Float32[40.0];
            comp_vmax = true,
        )
        @test all(iszero, photos.vmax)
        @test all(iszero, photos.agd)
        @test all(iszero, photos.rd)
        @test all(iszero, photos.adt)
        @test all(iszero, photos.adtmm)

        photos.tstress .= 1.0f0
        photos.lambda .= 0.8f0
        photos.vmax .= 1.0f0
        photosynthesis_C3!(
            cft1,
            photos,
            Float32[0.0],
            Float32[12.0],
            Float32[20.0],
            Float32[40.0];
            comp_vmax = false,
        )
        @test all(iszero, photos.adt)
        @test all(iszero, photos.adtmm)
    end

    @testset "C4 photosynthesis CPU smoke test" begin
        _, _, _, photos = init_crop(1, identity)
        photos.tstress .= 1.0f0

        photosynthesis_C4!(
            cft3,
            photos,
            Float32[10.0],
            Float32[12.0],
            Float32[25.0];
            comp_vmax = true,
        )

        @test photos.lambda == Float32[0.4]
        @test all(isfinite, photos.vmax)
        @test all(isfinite, photos.agd)
        @test all(isfinite, photos.rd)
        @test all(isfinite, photos.adt)
        @test all(isfinite, photos.adtmm)
        @test all(photos.vmax .>= 0.0f0)
        @test all(photos.agd .>= 0.0f0)
        @test all(photos.adt .>= 0.0f0)

        photos.tstress .= 0.0f0
        photosynthesis_C4!(
            cft3,
            photos,
            Float32[10.0],
            Float32[12.0],
            Float32[25.0];
            comp_vmax = true,
        )
        @test all(iszero, photos.vmax)
        @test all(iszero, photos.agd)
        @test all(iszero, photos.rd)
        @test all(iszero, photos.adt)
        @test all(iszero, photos.adtmm)
    end
end
