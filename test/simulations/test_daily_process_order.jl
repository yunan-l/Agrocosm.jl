using Test

@testset "C3/C4 daily process-order contract" begin
    driver = "daily_crop.jl"
    source = read(joinpath(@__DIR__, "..", "..", "src", "simulations", driver), String)
    daily_start = findfirst("function _daily_crop!", source)
    daily_start === nothing && error("missing common daily crop driver")
    source = source[first(daily_start):end]

    position = (call, start = firstindex(source)) -> begin
        pattern = Regex("(?m)^[ \\t]*" * call * "[ \\t]*\\(")
        location = findnext(pattern, source, start)
        location === nothing && error("missing $call call in $driver")
        first(location)
    end

    phenology_position = position("phenology_crop!")
    post_phenology_apar = position("_pathway_apar!", phenology_position)
    post_phenology_temp_stress = position("temp_stress", phenology_position)
    post_phenology_photosynthesis = position("photosynthesis!", phenology_position)

    @test position("update_climbuf!") < position("cultivate!")
    @test position("litter_tillage!") < position("tillage_hydraulics!") <
          position("pedotransfer!")
    @test position("_pathway_albedo!") < position("petpar!") < position("snow!")
    @test position("pedotransfer!") <
          position("update_surface_litter_properties!") <
          position("soil_temperature!")
    @test position("soil_temperature!") < position("soil_cn_decomposition!")
    @test position("soil_cn_decomposition!") < position("nitrogen_deposition!") <
          phenology_position
    @test position("_pathway_apar!") < position("temp_stress") <
          position("photosynthesis!") <
          position("prepare_prephenology_canopy_conductance!") < phenology_position
    @test position("cultivate!") < phenology_position < position("harvest_crop!")
    @test position("harvest_crop!") < position("route_harvest_residues!") <
          position("interception!") < position("add_snowmelt_to_precipitation!") <
          position("soil_infiltration!")
    # The rain/melt routing is common to both nitrogen-limit modes. Keeping
    # this call at loop scope prevents the legacy (`false`) path from silently
    # retaining the former pre-interception snowmelt behavior.
    @test occursin(
        r"(?m)^        add_snowmelt_to_precipitation!\($",
        source,
    )
    @test position("add_snowmelt_to_precipitation!") <
          position("record_water_balance_after_snow!") <
          position("soil_infiltration!")
    @test position("soil_infiltration!") < post_phenology_apar <
          post_phenology_temp_stress < post_phenology_photosynthesis <
          position("transpiration!") < position("solve_lambda!")
    @test position("solve_lambda!") < position("crop_carbon!") <
          position("terminate_failed_crop!") < position("evaporation!")
    @test position("solve_lambda!") < position("acquire_crop_nitrogen!") <
          position("limit_vcmax_by_nitrogen!") <
          position("recouple_nitrogen_water!") <
          position("finalize_nitrogen_limited_transpiration!") <
          position("crop_carbon!") <
          position("allocate_crop_nitrogen!")
    @test position("evaporation!") <
          position("soil_evapotranspiration!") < position("post_crop_nitrogen_losses!")
end
