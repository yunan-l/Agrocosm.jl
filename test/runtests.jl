using Agrocosm
using Test

@testset "Agrocosm.jl" begin
    include("numerics/test_lpj_bisect.jl")
    include("processes/initialization/test_initialization.jl")
    include("processes/climate/test_snow.jl")
    include("processes/crop/test_photosynthesis.jl")
    include("processes/crop/test_lambda_solver_c3.jl")
    include("processes/crop/test_lambda_solver_c4.jl")
    include("processes/crop/test_lambda_water_coupling.jl")
    include("processes/soil/test_soil_water.jl")
    include("diagnostics/test_water_balance.jl")
end
