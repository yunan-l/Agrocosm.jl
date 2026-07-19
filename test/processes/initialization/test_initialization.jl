using Agrocosm
using Test

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
