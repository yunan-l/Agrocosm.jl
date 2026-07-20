using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

# Defines the shared deterministic CPU/GPU precision fixture and also runs its
# CPU-only regression before the GPU comparisons below.
include("test_daily_crop_C3_precision.jl")

@testset "CUDA C3 Float32/Float64 precision support" begin
    for T in (Float32, Float64)
        cpu = run_c3_precision_smoke(T, identity)
        gpu = run_c3_precision_smoke(T, CuArray)

        gpu_npp = Array(gpu.output.crop.npp)
        gpu_water = Array(gpu.soil.water.storage)
        @test eltype(gpu.output.crop.npp) == T
        @test eltype(gpu.soil.water.storage) == T
        @test all(isfinite, gpu_npp)
        @test all(isfinite, gpu_water)
        @test gpu_npp ≈ cpu.output.crop.npp rtol = T(5e-5) atol = T(5e-6)
        @test gpu_water ≈ cpu.soil.water.storage rtol = T(5e-5) atol = T(5e-6)
    end
end
