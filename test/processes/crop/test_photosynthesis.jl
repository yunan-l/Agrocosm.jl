using Agrocosm
using Test

@testset "C3 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    photos = crop.photosynthesis
    photos.temperature_stress .= 1.0f0

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
    @test all(isfinite, photos.gross_assimilation)
    @test all(isfinite, photos.leaf_respiration)
    @test all(isfinite, photos.net_assimilation)
    @test all(isfinite, photos.water_limited_assimilation)
    @test all(photos.vmax .>= 0.0f0)
    @test all(photos.gross_assimilation .>= 0.0f0)
    @test all(photos.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
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
    @test all(iszero, photos.gross_assimilation)
    @test all(iszero, photos.leaf_respiration)
    @test all(iszero, photos.net_assimilation)
    @test all(iszero, photos.water_limited_assimilation)

    photos.temperature_stress .= 1.0f0
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
    @test all(iszero, photos.net_assimilation)
    @test all(iszero, photos.water_limited_assimilation)
end

@testset "C4 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    photos = crop.photosynthesis
    photos.temperature_stress .= 1.0f0

    photosynthesis_C4!(
        cft3,
        photos,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vmax = true,
    )

    @test photos.lambda == Float32[0.8]
    @test all(isfinite, photos.vmax)
    @test all(isfinite, photos.gross_assimilation)
    @test all(isfinite, photos.leaf_respiration)
    @test all(isfinite, photos.net_assimilation)
    @test all(isfinite, photos.water_limited_assimilation)
    @test all(photos.vmax .>= 0.0f0)
    @test all(photos.gross_assimilation .>= 0.0f0)
    @test all(photos.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
    photosynthesis_C4!(
        cft3,
        photos,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vmax = true,
    )
    @test all(iszero, photos.vmax)
    @test all(iszero, photos.gross_assimilation)
    @test all(iszero, photos.leaf_respiration)
    @test all(iszero, photos.net_assimilation)
    @test all(iszero, photos.water_limited_assimilation)
end
