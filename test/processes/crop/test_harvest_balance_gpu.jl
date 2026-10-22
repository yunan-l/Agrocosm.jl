using CUDA

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

include("test_harvest_balance.jl")

@testset "CUDA last-harvest annual yield" begin
    test_last_harvest_yield(CuArray)
end
