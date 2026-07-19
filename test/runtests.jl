using Agrocosm
using Test

@testset "Agrocosm.jl" begin
    include("numerics/test_lpj_bisect.jl")
    include("processes/initialization/test_initialization.jl")
    include("processes/climate/test_snow.jl")
    include("processes/climate/test_readclimate.jl")
    include("processes/crop/test_photosynthesis.jl")
    include("processes/crop/test_lambda_solver_c3.jl")
    include("processes/crop/test_lambda_solver_c4.jl")
    include("processes/crop/test_lambda_water_coupling.jl")
    include("processes/crop/test_fertilizer.jl")
    include("processes/crop/test_nitrogen_allocation.jl")
    include("processes/crop/test_nitrogen_demand.jl")
    include("processes/crop/test_nitrogen_uptake.jl")
    include("processes/soil/test_mineralization_immobilization.jl")
    include("processes/soil/test_nitrification.jl")
    include("processes/soil/test_denitrification.jl")
    include("processes/soil/test_volatilization.jl")
    include("processes/soil/test_soil_water.jl")
    include("diagnostics/test_water_balance.jl")
    include("diagnostics/test_nitrogen_balance.jl")
end
