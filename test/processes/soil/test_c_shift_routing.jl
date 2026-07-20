using Agrocosm
using Test

@testset "Fixed post-spin-up c_shift routing" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    fast_shift = Float32[0.40, 0.25, 0.15, 0.10, 0.10]
    slow_shift = Float32[0.55, 0.20, 0.10, 0.10, 0.05]

    soil.carbon.shift_fast[:, 1] .= fast_shift
    soil.carbon.shift_slow[:, 1] .= slow_shift
    soil.nitrogen.shift_fast[:, 1] .= fast_shift
    soil.nitrogen.shift_slow[:, 1] .= slow_shift
    soil.carbon.litter[Agrocosm.ROOT_LITTER, 1] = 10.0f0
    soil.nitrogen.litter[Agrocosm.ROOT_LITTER, 1] = 0.4f0
    soil.carbon.litter_response .= 0.0f0
    soil.nitrogen.litter_response .= 0.0f0
    soil.carbon.litter_response[Agrocosm.ROOT_LITTER] = 1.0f0
    soil.nitrogen.litter_response[Agrocosm.ROOT_LITTER] = 1.0f0
    soil.thermal.temperature .= 10.0f0
    soil.water.relative_content .= 0.5f0
    soil.water.holding_capacity_storage .= 80.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.nitrogen.ammonium .= 0.2f0
    soil.nitrogen.nitrate .= 0.1f0

    @test sum(soil.carbon.shift_fast; dims = 1)[1] ≈ 1.0f0
    @test sum(soil.carbon.shift_slow; dims = 1)[1] ≈ 1.0f0

    soil_carbon!(crop, soil)
    soil_nitrogen!(crop, soil; air_temperature = Float32[10], wind_speed = Float32[1.5])

    decomposed_carbon = sum(soil.carbon.decomposed_litter)
    decomposed_nitrogen = sum(soil.nitrogen.decomposed_litter)
    retained_carbon = decomposed_carbon * (1.0f0 - lpjmlparams.atmfrac)
    retained_nitrogen = decomposed_nitrogen * (1.0f0 - lpjmlparams.atmfrac)

    @test soil.carbon.litter_to_fast[:, 1] ≈
          fast_shift .* retained_carbon .* lpjmlparams.fastfrac
    @test soil.carbon.litter_to_slow[:, 1] ≈
          slow_shift .* retained_carbon .* (1.0f0 - lpjmlparams.fastfrac)
    @test soil.nitrogen.litter_to_fast[:, 1] ≈
          fast_shift .* retained_nitrogen .* lpjmlparams.fastfrac
    @test soil.nitrogen.litter_to_slow[:, 1] ≈
          slow_shift .* retained_nitrogen .* (1.0f0 - lpjmlparams.fastfrac)
    @test sum(soil.carbon.litter_to_fast .+ soil.carbon.litter_to_slow) ≈
          retained_carbon atol = 1.0f-6
    @test sum(soil.nitrogen.litter_to_fast .+ soil.nitrogen.litter_to_slow) ≈
          retained_nitrogen atol = 1.0f-7
    @test sum(root_distribution(0.96f0)) ≈ 1.0 atol = 1.0e-6
end
