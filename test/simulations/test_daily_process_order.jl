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
    @test position("cultivate!") < position("phenology_crop!") < position("harvest_crop!")
    @test position("harvest_crop!") < position("route_harvest_residues!")
    @test position("crop_carbon!") < position("terminate_failed_crop!")
    @test position("terminate_failed_crop!") < position("post_crop_nitrogen_losses!")
end
