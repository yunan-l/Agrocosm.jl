using Agrocosm
using Test

@testset "Agrocosm.jl" begin
    include("test_initialization.jl")
    include("test_photosynthesis.jl")
    include("test_water_balance.jl")
end
