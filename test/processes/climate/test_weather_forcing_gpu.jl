using Agrocosm
using CUDA
using Test

CUDA.functional() || error("This targeted CUDA check requires a working NVIDIA device")
CUDA.allowscalar(false)
@testset "Weather forcing CPU/CUDA equivalence" begin
    forcing = reshape(Float32.(1:40), 4, 2, 5)
    cpu = Agrocosm.init_weather(Float32, 2, identity)
    gpu = Agrocosm.init_weather(Float32, 2, CUDA.cu)
    Agrocosm.apply_weather_forcing!(cpu, forcing, 3)
    Agrocosm.apply_weather_forcing!(gpu, CUDA.cu(forcing), 3)
    CUDA.synchronize()
    for field in (:temp, :prec, :swr, :lwr, :wind, :no3_deposition, :nh4_deposition)
        @test Array(getproperty(gpu, field)) == getproperty(cpu, field)
    end
end
