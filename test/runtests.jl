using Agrocosm
using Test

# During the Terrarium migration the crop/soil physics is being re-expressed as
# continuous-time Terrarium processes (plan Phases 2–6). The legacy test files
# under test/processes/**, test/diagnostics/**, and test/simulations/** target
# the deleted standalone API and are re-enabled as each process is ported. Only
# the infrastructure-free tests run in the interim.
@testset "Agrocosm.jl" begin
    @testset "Numerics" begin
        include("numerics/test_lpj_bisect.jl")
    end

    @testset "Crop parameters" begin
        include("processes/crop/test_pft_registry.jl")
    end

    @testset "Crop processes" begin
        include("crop/test_root_distribution.jl")
        include("crop/test_photosynthesis.jl")
        include("crop/test_stomatal_conductance.jl")
        include("crop/test_carbon_dynamics.jl")
        include("crop/test_phenology.jl")
        include("crop/test_nitrogen_limitation.jl")
        include("crop/test_nitrogen_demand.jl")
        include("crop/test_plant_available_water.jl")
        include("crop/test_nitrogen_uptake.jl")
        include("crop/test_soil_decomposition_response.jl")
        include("crop/test_growth_respiration.jl")
        include("crop/test_harvest_index.jl")
    end
end
