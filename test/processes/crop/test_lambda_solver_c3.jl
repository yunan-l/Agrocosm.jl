using Agrocosm
using Test

@testset "Backend-compatible C3 lambda kernel" begin
    crop, _, _, photos = init_crop(2, identity)
    pet = init_pet(2, identity)

    target_lambda = 0.5f0
    vmax = 2.0f0
    tstress = 1.0f0
    co2 = Float32[40.0, 40.0]
    temp = Float32[20.0, 20.0]
    apar = 1.0f6
    daylength = 12.0f0
    fpar = 0.8f0

    target_adtmm = Agrocosm.c3_adtmm_scalar(
        target_lambda, vmax, tstress, cft1.b, co2[1], temp[1], apar, daylength,
    )
    fac = target_adtmm / (1.0f0 - target_lambda)
    gpd = fac * 1.6f0 / (co2[1] * 1.0f-5)
    target_conductance = gpd / (daylength * 3600.0f0) + cft1.gmin * fpar

    photos.vmax .= vmax
    photos.tstress .= tstress
    crop.apar .= apar
    crop.fpar .= fpar
    crop.gp .= Float32[target_conductance, 0.0]
    pet.daylength .= daylength

    solve_lambda_c3!(cft1, photos, crop, pet, temp, co2)

    @test photos.lambda[1] ≈ target_lambda atol = 2.0f-3
    @test photos.lambda[2] == 0.0f0
end

@testset "LPJ-compatible C3 lambda solver" begin
    vmax = 2.0f0
    tstress = 1.0f0
    b = cft1.b
    co2 = 40.0f0
    temp = 20.0f0
    apar = 1.0f6
    daylength = 12.0f0
    target_lambda = 0.5f0

    target_adtmm = Agrocosm.c3_adtmm_scalar(
        target_lambda, vmax, tstress, b, co2, temp, apar, daylength,
    )
    @test target_adtmm > 0.0f0

    _, _, _, photos = init_crop(1, identity)
    photos.tstress .= tstress
    photos.lambda .= target_lambda
    photos.vmax .= vmax
    photosynthesis_C3!(
        cft1,
        photos,
        Float32[apar],
        Float32[daylength],
        Float32[temp],
        Float32[co2];
        comp_vmax = false,
    )
    @test target_adtmm ≈ photos.adtmm[1] atol = 1.0f-6

    fac = target_adtmm / (1.0f0 - target_lambda)

    lambda, iterations, residual = solve_lambda_c3_lpj(
        fac, vmax, tstress, b, co2, temp, apar, daylength,
    )

    @test 0.02f0 <= lambda <= 0.85f0
    @test lambda ≈ target_lambda atol = 2.0f-3
    @test abs(residual) < 1.0f-3
    @test iterations <= 30
end
