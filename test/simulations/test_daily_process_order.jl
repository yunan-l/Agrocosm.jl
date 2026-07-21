using Test

@testset "C3/C4 daily process-order contract" begin
    for driver in ("daily_crop_C3.jl", "daily_crop_C4.jl")
        source = read(joinpath(@__DIR__, "..", "..", "src", "simulations", driver), String)

        position = call -> begin
            match = findfirst("\n        " * call * "(", source)
            match === nothing && error("missing $call call in $driver")
            first(match)
        end

        @test position("update_climbuf!") < position("cultivate!")
        @test position("albedo!") < position("petpar!") < position("snow!")
        @test position("pedotransfer!") <
              position("update_surface_litter_properties!") <
              position("soil_temperature!")
    end
end
