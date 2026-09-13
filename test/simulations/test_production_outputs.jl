using Agrocosm
using Test

@testset "Every production output variable has a runtime contract" begin
    # `production_output_variables()` is built in the global entry script, and
    # `OutputVariable` refuses a field with no entry in the runtime-contract
    # registry. Nothing in the unit suite called it, so a variable added to the
    # writer without its contract passed every local test and then aborted 32
    # MPI ranks on the first real job. This closes that path: the list is
    # constructed here, which is the same call the entry makes.
    include(joinpath(@__DIR__, "..", "..", "scripts", "run_global_wheat_cpu.jl"))
    variables = production_output_variables()
    @test !isempty(variables)
    for variable in variables
        # `output_variable_spec` returns (spec, frequency).
        spec, frequency = Agrocosm.output_variable_spec(variable.group, variable.field)
        @test !isempty(spec.units)
        @test !isempty(spec.description)
        @test frequency in (:annual, :daily)
    end
    # The one this test was written for, and the reason it is annual: a season
    # accumulation emitted at the calendar-year boundary, not a daily value.
    lodging = only(v for v in variables if v.field === :lodging_exposure)
    @test lodging.group === :crop
    @test last(Agrocosm.output_variable_spec(:crop, :lodging_exposure)) === :annual
end
