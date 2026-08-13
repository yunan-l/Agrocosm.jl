using Agrocosm
using Enzyme
using Test

@testset "Agrocosm Enzyme adapter" begin
    include("ad/test_adapter_contract.jl")
    include("ad/test_enzyme_adapter.jl")
    include("ad/test_enzyme_daily_transition.jl")
    include("ad/test_enzyme_seasonal_loss.jl")
    include("ad/test_enzyme_reverse.jl")
    include("ad/test_enzyme_model_parameter_reverse.jl")
    include("ad/test_enzyme_365day.jl")
    include("ad/test_enzyme_blockwise.jl")
end
