using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA C4 lambda solver" begin
    crop, _, _, photos = init_crop(2, CuArray)
    pet = init_pet(2, CuArray)

    target_lambda = 0.2f0
    vmax = 2.0f0
    tstress = 1.0f0
    co2_cpu = Float32[40.0, 40.0]
    temp_cpu = Float32[25.0, 25.0]
    apar = 1.0f6
    daylength = 12.0f0
    fpar = 0.8f0

    target_adtmm = Agrocosm.c4_adtmm_scalar(
        target_lambda, vmax, tstress, cft3.b, temp_cpu[1], apar, daylength,
    )
    fac = target_adtmm / (1.0f0 - target_lambda)
    gpd = fac * 1.6f0 / (co2_cpu[1] * 1.0f-5)
    target_conductance = gpd / (daylength * 3600.0f0) + cft3.gmin * fpar

    photos.vmax .= vmax
    photos.tstress .= tstress
    crop.apar .= apar
    crop.fpar .= fpar
    crop.gp .= CuArray(Float32[target_conductance, 0.0])
    pet.daylength .= daylength

    solve_lambda_c4!(
        cft3, photos, crop, pet, CuArray(temp_cpu), CuArray(co2_cpu),
    )

    lambda_cpu = Array(photos.lambda)
    @test lambda_cpu[1] ≈ target_lambda atol = 2.0f-3
    @test lambda_cpu[2] == 0.0f0
end
