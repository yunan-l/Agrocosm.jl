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

    for field in (:sowing_callback, :harvest_callback, :harvesting_mask)
        @test Array(getproperty(gpu.output.calendar, field)) ==
            Array(getproperty(cpu.output.calendar, field))
    end
    @test callback_days(gpu.output.calendar.sowing_callback) == [100, 465]
    @test callback_days(gpu.output.calendar.harvest_callback) == [137, 502]
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
    for field in (:harvesting_mask, :sowing_callback, :harvest_callback)
        @test all(
            iszero,
            Array(getproperty(gpu.output.calendar, field))[inactive_days, :],
        )
    end

    for (container, fields) in (
        (gpu.crop.phenology, (:vdsum, :husum, :fphu, :growing_days, :is_growing)),
        (gpu.crop.canopy,
         (:lai, :flaimax, :laimax_adjusted, :lai_npp_deficit,
          :phenology_fraction, :fpar, :apar)),
        (gpu.crop.carbon,
         (:biomass, :leaf, :root, :pool, :storage, :organs, :yield,
          :npp, :respiration)),
        (gpu.crop.nitrogen,
         (:total, :uptake, :auto_fertilizer, :leaf, :root, :pool, :storage,
          :demand_total, :demand_leaf, :pending_manure, :pending_fertilizer,
          :seed_input, :prescribed_manure_input,
          :prescribed_fertilizer_input, :harvest_export, :stress_sum,
          :stress, :deficit)),
        (gpu.crop.water,
         (:canopy_conductance, :transpiration, :canopy_wet, :interception,
          :transpiration_layer, :deficit, :demand_sum, :supply_sum, :stress,
          :waterlogging_days)),
        (gpu.crop.photosynthesis,
         (:gross_assimilation, :net_assimilation,
          :water_limited_assimilation, :leaf_respiration, :potential_vmax,
          :vmax, :nitrogen_limitation, :lambda)),
    )
        for field in fields
            @test all(iszero, Array(getproperty(container, field)))
        end
    end
    @test Array(gpu.crop.phenology.harvesting) == Bool[true]
    @test Array(gpu.crop.phenology.harvesting_previous) == Bool[true]
    @test Array(gpu.crop.water.waterlogging_stress) == Float32[1]
end
