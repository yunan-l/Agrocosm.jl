using Agrocosm
using Test

@testset "C3 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    photos = crop.auxiliary.photosynthesis
    photos.temperature_stress .= 1.0f0

    photosynthesis_C3!(
        cft1,
        crop,
        Float32[10.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = true,
    )

    @test photos.lambda == Float32[0.8]
    @test all(isfinite, photos.vcmax)
    @test all(isfinite, crop.fluxes.carbon.gross_assimilation)
    @test all(isfinite, crop.fluxes.carbon.leaf_respiration)
    @test all(isfinite, crop.fluxes.carbon.net_assimilation)
    @test all(isfinite, crop.fluxes.carbon.water_limited_assimilation)
    @test all(photos.vcmax .>= 0.0f0)
    @test all(crop.fluxes.carbon.gross_assimilation .>= 0.0f0)
    @test all(crop.fluxes.carbon.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
    photosynthesis_C3!(
        cft1,
        crop,
        Float32[10.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = true,
    )
    @test all(iszero, photos.vcmax)
    @test all(iszero, crop.fluxes.carbon.gross_assimilation)
    @test all(iszero, crop.fluxes.carbon.leaf_respiration)
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)

    photos.temperature_stress .= 1.0f0
    photos.lambda .= 0.8f0
    photos.vcmax .= 1.0f0
    photosynthesis_C3!(
        cft1,
        crop,
        Float32[0.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = false,
    )
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)
end

@testset "C4 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    photos = crop.auxiliary.photosynthesis
    photos.temperature_stress .= 1.0f0

    photosynthesis_C4!(
        cft3,
        crop,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vcmax = true,
    )

    @test photos.lambda == Float32[0.8]
    @test all(isfinite, photos.vcmax)
    @test all(isfinite, crop.fluxes.carbon.gross_assimilation)
    @test all(isfinite, crop.fluxes.carbon.leaf_respiration)
    @test all(isfinite, crop.fluxes.carbon.net_assimilation)
    @test all(isfinite, crop.fluxes.carbon.water_limited_assimilation)
    @test all(photos.vcmax .>= 0.0f0)
    @test all(crop.fluxes.carbon.gross_assimilation .>= 0.0f0)
    @test all(crop.fluxes.carbon.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
    photosynthesis_C4!(
        cft3,
        crop,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vcmax = true,
    )
    @test all(iszero, photos.vcmax)
    @test all(iszero, crop.fluxes.carbon.gross_assimilation)
    @test all(iszero, crop.fluxes.carbon.leaf_respiration)
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)
end
