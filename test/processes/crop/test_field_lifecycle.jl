using Agrocosm
using Test

include("../../helpers/crop_lifecycle_fixture.jl")

field_arrays(container) =
    [getproperty(container, field) for field in fieldnames(typeof(container))]

function daily_owned_arrays(crop)
    return vcat(
        field_arrays(crop.fluxes.carbon),
        field_arrays(crop.fluxes.nitrogen),
        field_arrays(crop.fluxes.water),
        [
            crop.auxiliary.phenology.fphu,
            crop.auxiliary.canopy.flaimax,
            crop.auxiliary.canopy.actual_lai,
            crop.auxiliary.canopy.albedo,
            crop.auxiliary.canopy.fpar,
            crop.auxiliary.canopy.apar,
            crop.auxiliary.canopy.canopy_conductance,
            crop.auxiliary.canopy.canopy_wet,
        ],
        field_arrays(crop.auxiliary.photosynthesis),
        [
            crop.auxiliary.stress.nitrogen_demand_total,
            crop.auxiliary.stress.nitrogen_demand_leaf,
            crop.auxiliary.stress.nitrogen_deficit,
            crop.auxiliary.stress.water_deficit,
            crop.auxiliary.root.zone_available_water,
        ],
        field_arrays(crop.workspace),
    )
end

function poison_daily_fields!(crop, value)
    crop.events.sowing .= 1
    crop.events.harvest .= 1
    for values in daily_owned_arrays(crop)
        values .= value
    end
    return nothing
end

function compare_nested_arrays(left, right)
    @test typeof(left) == typeof(right)
    for field in fieldnames(typeof(left))
        left_value = getfield(left, field)
        right_value = getfield(right, field)
        if left_value isa AbstractArray
            @test left_value == right_value
        elseif isstructtype(typeof(left_value))
            compare_nested_arrays(left_value, right_value)
        else
            @test isequal(left_value, right_value)
        end
    end
    return nothing
end

function run_owner_overwrite_case(::Type{T}; poison::Bool) where {T <: AbstractFloat}
    initial_data = lifecycle_initial_data(T)
    initial_data.ModelState.crop.sdate .= Int32(1)
    simulation = initialize_simulation(
        cft1, initial_data;
        T = T,
        days = 2,
        diagnostics = false,
        auto_fertilizer = false,
    )
    one_day = lifecycle_climate(T, 1)
    run_simulation!(simulation, one_day; spinup = false)
    poison && poison_daily_fields!(simulation.crop, T(123))
    run_simulation!(simulation, one_day; spinup = false)
    return simulation
end

@testset "Owner processes overwrite daily fields without a global reset" begin
    for T in (Float32, Float64)
        clean = run_owner_overwrite_case(T; poison = false)
        poisoned = run_owner_overwrite_case(T; poison = true)

        compare_nested_arrays(clean.crop, poisoned.crop)
        compare_nested_arrays(clean.soil, poisoned.soil)
        compare_nested_arrays(clean.output, poisoned.output)
    end
end

@testset "State and auxiliary fields have an explicit cross-day contract" begin
    crop = init_crop(2, identity)
    @test propertynames(crop.auxiliary.stress) == (
        :nitrogen_demand_total,
        :nitrogen_demand_leaf,
        :nitrogen_deficit,
        :water_deficit,
    )
    @test crop.auxiliary.root.distribution isa AbstractVector
    @test :canopy ∈ propertynames(crop.state)
    @test :lai ∈ propertynames(crop.state.canopy)
    @test :flaimax ∈ propertynames(crop.auxiliary.canopy)
    @test :phenology_fraction ∉ propertynames(crop.state.canopy)
    @test :actual_lai ∈ propertynames(crop.auxiliary.canopy)
    @test propertynames(crop.state.phenology) == (
        :vdsum, :husum, :senescence, :senescence_previous,
        :harvesting, :harvesting_previous, :growing_days, :is_growing,
    )
    @test propertynames(crop.auxiliary.phenology) == (:phu, :winter_type, :fphu)
    @test propertynames(crop.auxiliary.calendar) == (:sowing_date,)
    @test :sufficiency ∈ propertynames(crop.state.nitrogen)
    @test :sufficiency ∈ propertynames(crop.state.water)
end

@testset "Workspace is outside scientific output and restart" begin
    crop = init_crop(1, identity)
    output = init_output(1, identity)
    restart = crop_restart_payload(crop)

    @test propertynames(restart) == (:state, :process_memory)
    @test restart.state === crop.state
    @test restart.state.canopy.lai_npp_deficit === crop.state.canopy.lai_npp_deficit
    @test :fphu ∉ propertynames(restart.state.phenology)
    @test restart.process_memory.calendar.sowing_date ===
          crop.auxiliary.calendar.sowing_date
    @test :workspace ∉ propertynames(restart)
    @test :workspace ∉ fieldnames(typeof(output))
    @test propertynames(output.annual) == (:yield, :harvest_date)
    @test isempty(fieldnames(typeof(crop.workspace)))
end
