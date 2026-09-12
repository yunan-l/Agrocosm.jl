using Agrocosm
using Enzyme
using Test

run_full_suite = "--full" in ARGS

# `test_zero_leaf_carbon_gradient.jl` builds a crop state directly, so it needs
# the same fixture helpers the unit suite uses.
include("helpers/model_state_fixture.jl")

@testset "Agrocosm Enzyme adapter" begin
    include("ad/test_zero_leaf_carbon_gradient.jl")
    include("ad/test_cft_offset_gradient.jl")
    include("ad/test_process_parameter_gradient.jl")
    include("ad/test_absolute_threshold_gradient.jl")
    include("ad/test_enzyme_daily_transition.jl")
    include("ad/test_enzyme_seasonal_loss.jl")
    include("ad/test_enzyme_365day.jl")
    if run_full_suite
        include("ad/test_adapter_contract.jl")
        include("ad/test_enzyme_adapter.jl")
        include("ad/test_enzyme_nitrogen_limit.jl")
        include("ad/test_enzyme_irrigation.jl")
        include("ad/test_enzyme_management_adaptation.jl")
        include("ad/test_enzyme_reverse.jl")
        include("ad/test_enzyme_model_parameter_reverse.jl")
        include("ad/test_enzyme_blockwise.jl")
    end
end
