using Test

@testset "C3/C4 daily process-order contract" begin
    driver = "daily_crop.jl"
    source = read(joinpath(@__DIR__, "..", "..", "src", "simulations", driver), String)
    daily_start = findfirst("function _daily_crop!", source)
    daily_start === nothing && error("missing common daily crop driver")
    source = source[first(daily_start):end]

    position = call -> begin
        pattern = Regex("(?m)^[ \\t]*" * call * "[ \\t]*\\(")
        location = findfirst(pattern, source)
        location === nothing && error("missing $call call in $driver")
        first(location)
    end

    @test position("update_climbuf!") < position("cultivate!")
    @test position("litter_tillage!") < position("tillage_hydraulics!") <
          position("pedotransfer!")
    @test position("_pathway_albedo!") < position("petpar!") < position("snow!")
    @test position("pedotransfer!") <
          position("update_surface_litter_properties!") <
          position("soil_temperature!")
end
