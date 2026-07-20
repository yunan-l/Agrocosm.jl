using Agrocosm
using Test

include("../helpers/crop_lifecycle_fixture.jl")

function lifecycle_simulation(days)
    return initialize_simulation(
        cft1,
        lifecycle_initial_data(Float32);
        T = Float32,
        days = days,
        diagnostics = false,
        auto_fertilizer = false,
    )
end

@testset "Annual climate blocks preserve crop events and soil memory" begin
    days = 730
    continuous = lifecycle_simulation(days)
    chunked = lifecycle_simulation(days)
    full_climate = lifecycle_climate(Float32, days)
    first_year = NamedTuple(
        field => getproperty(full_climate, field) isa AbstractMatrix ?
            getproperty(full_climate, field)[1:365, :] : getproperty(full_climate, field)
        for field in propertynames(full_climate)
    )
    second_year = NamedTuple(
        field => getproperty(full_climate, field) isa AbstractMatrix ?
            getproperty(full_climate, field)[366:730, :] : getproperty(full_climate, field)
        for field in propertynames(full_climate)
    )

    run_simulation!(continuous, full_climate; spinup = false)
    run_simulation!(chunked, [first_year, second_year]; spinup = false)

    @test chunked.simulated_days == days
    for field in (:sowing_callback, :harvest_callback, :harvesting_mask)
        @test getproperty(chunked.output.calendar, field) ==
            getproperty(continuous.output.calendar, field)
    end
    for field in (:fphu, :biomass, :lai, :npp)
        @test getproperty(chunked.output.crop, field) ≈
            getproperty(continuous.output.crop, field)
    end
    for (chunked_state, continuous_state) in (
        (chunked.soil.water.storage, continuous.soil.water.storage),
        (chunked.soil.water.ice_storage, continuous.soil.water.ice_storage),
        (chunked.soil.thermal.temperature, continuous.soil.thermal.temperature),
        (chunked.soil.carbon.litter, continuous.soil.carbon.litter),
        (chunked.soil.carbon.fast, continuous.soil.carbon.fast),
        (chunked.soil.carbon.slow, continuous.soil.carbon.slow),
        (chunked.soil.nitrogen.litter, continuous.soil.nitrogen.litter),
        (chunked.soil.nitrogen.fast, continuous.soil.nitrogen.fast),
        (chunked.soil.nitrogen.slow, continuous.soil.nitrogen.slow),
        (chunked.soil.nitrogen.nitrate, continuous.soil.nitrogen.nitrate),
        (chunked.soil.nitrogen.ammonium, continuous.soil.nitrogen.ammonium),
    )
        @test chunked_state ≈ continuous_state
    end
    @test callback_days(chunked.output.calendar.sowing_callback) == [100, 465]
    @test callback_days(chunked.output.calendar.harvest_callback) == [137, 502]
end

