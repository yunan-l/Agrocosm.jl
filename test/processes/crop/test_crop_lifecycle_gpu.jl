using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

include("../../helpers/crop_lifecycle_fixture.jl")

@testset "CUDA two-year prescribed crop lifecycle" begin
    cpu = run_lifecycle_fixture(identity)
    gpu = run_lifecycle_fixture(CuArray)
    synchronize()

    for field in (:sowing_event, :harvest_event, :harvesting_mask)
        @test Array(getproperty(gpu.output.calendar, field)) ==
            Array(getproperty(cpu.output.calendar, field))
    end
    @test event_days(gpu.output.calendar.sowing_event) == [100, 465]
    @test event_days(gpu.output.calendar.harvest_event) == [137, 502]
    @test Array(gpu.output.crop.growing_mask) == Array(cpu.output.crop.growing_mask)
    @test Array(gpu.output.crop.fphu) ≈ Array(cpu.output.crop.fphu) rtol = 2.0f-5

    inactive_days = vcat(138:464, 503:730)
    for field in (
        :gpp, :npp, :lambda, :potential_vmax, :vmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu,
    )
        @test all(iszero, Array(getproperty(gpu.output.crop, field))[inactive_days, :])
    end
    @test all(iszero, Array(gpu.output.crop.growing_mask)[inactive_days, :])
    for field in (:harvesting_mask, :sowing_event, :harvest_event)
        @test all(
            iszero,
            Array(getproperty(gpu.output.calendar, field))[inactive_days, :],
        )
    end

    for (container, fields) in (
        (gpu.crop.state.phenology, (:vdsum, :husum, :fphu, :growing_days, :is_growing)),
        (gpu.crop.state.canopy,
         (:lai, :laimax_adjusted, :lai_npp_deficit)),
        (gpu.crop.auxiliary.canopy,
         (:flaimax, :phenology_fraction, :fpar, :apar,
          :canopy_conductance, :canopy_wet)),
        (gpu.crop.state.carbon,
         (:biomass, :leaf, :root, :pool, :storage)),
        (gpu.crop.fluxes.carbon,
         (:yield, :npp, :respiration, :gross_assimilation, :net_assimilation,
          :water_limited_assimilation, :leaf_respiration)),
        (gpu.crop.state.nitrogen,
         (:total, :leaf, :root, :pool, :storage, :pending_manure,
          :pending_fertilizer, :stress_sum)),
        (gpu.crop.fluxes.nitrogen,
         (:uptake, :auto_fertilizer, :seed_input, :prescribed_manure_input,
          :prescribed_fertilizer_input, :harvest_export)),
        (gpu.crop.state.water,
         (:demand_sum, :supply_sum, :waterlogging_days)),
        (gpu.crop.fluxes.water,
         (:transpiration, :interception, :transpiration_layer)),
        (gpu.crop.auxiliary.stress,
         (:nitrogen_demand_total, :nitrogen_demand_leaf, :nitrogen,
          :nitrogen_deficit, :water_deficit, :water)),
        (gpu.crop.auxiliary.photosynthesis,
         (:potential_vmax, :vmax, :nitrogen_limitation, :lambda)),
    )
        for field in fields
            @test all(iszero, Array(getproperty(container, field)))
        end
    end
    @test Array(gpu.crop.state.phenology.harvesting) == Bool[true]
    @test Array(gpu.crop.state.phenology.harvesting_previous) == Bool[true]
    @test Array(gpu.crop.auxiliary.stress.waterlogging) == Float32[1]
end
