using Agrocosm
using Test

field_arrays(container) =
    [getproperty(container, field) for field in fieldnames(typeof(container))]

function daily_owned_arrays(state::ModelState)
    return vcat(
        field_arrays(state.fluxes.crop.carbon),
        field_arrays(state.fluxes.crop.nitrogen),
        field_arrays(state.fluxes.crop.water),
        [
            state.auxiliary.crop.phenology.fphu,
            state.auxiliary.crop.canopy.flaimax,
            state.auxiliary.crop.canopy.actual_lai,
            state.auxiliary.crop.canopy.albedo,
            state.auxiliary.crop.canopy.fpar,
            state.auxiliary.crop.canopy.apar,
            state.auxiliary.crop.canopy.canopy_conductance,
            state.auxiliary.crop.canopy.canopy_wet,
        ],
        # Every photosynthesis auxiliary EXCEPT `pmodel_chi`, which is written
        # only when the P-model is switched on. It is not daily-owned in the
        # sense this test means: with `pmodel_beta = 0` the kernel that writes it
        # never runs, and the solver reads its zero as "no bound" and falls back
        # to LPJmL's own 0.85. Poisoning it would therefore assert that an
        # inactive mechanism still writes every day, which is the opposite of the
        # ablation contract everything here ships with.
        [getfield(state.auxiliary.crop.photosynthesis, name)
         for name in fieldnames(typeof(state.auxiliary.crop.photosynthesis))
         if name !== :pmodel_chi],
        [
            state.auxiliary.crop.stress.nitrogen_demand_total,
            state.auxiliary.crop.stress.nitrogen_demand_leaf,
            state.auxiliary.crop.stress.nitrogen_deficit,
            state.auxiliary.crop.stress.water_deficit,
            state.auxiliary.crop.root.zone_available_water,
        ],
        field_arrays(state.workspace.crop),
    )
end

function poison_daily_fields!(state::ModelState, value)
    state.events.crop.sowing .= 1
    state.events.crop.harvest .= 1
    for values in daily_owned_arrays(state)
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
        fertilizer = :yes,
    )
    one_day = lifecycle_climate(T, 1)
    run_simulation!(simulation, one_day; spinup = false)
    poison && poison_daily_fields!(simulation.state, T(123))
    run_simulation!(simulation, one_day; spinup = false)
    return simulation
end

@testset "Owner processes overwrite daily fields without a global reset" begin
    for T in (Float32, Float64)
        clean = run_owner_overwrite_case(T; poison = false)
        poisoned = run_owner_overwrite_case(T; poison = true)

        compare_nested_arrays(clean.state, poisoned.state)
        compare_nested_arrays(clean.output, poisoned.output)
    end
end

@testset "State and auxiliary fields have an explicit cross-day contract" begin
    crop = init_crop(2, identity)
    @test propertynames(crop.auxiliary.photosynthesis) == (
        :potential_vcmax, :vcmax, :nitrogen_limitation, :lambda,
        :temperature_stress,
        # The P-model's least-cost optimal ci/ca, written only when the P-model
        # is on and read as the UPPER BOUND of the lambda solve. Not a
        # replacement for `lambda`: LPJmL's solve encodes soil water limitation
        # and the optimum encodes atmospheric demand, and measured on real cells
        # the two correlate 0.03 to 0.05, so neither stands in for the other.
        :pmodel_chi,
    )
    @test propertynames(crop.auxiliary.stress) == (
        :nitrogen_demand_total,
        :nitrogen_demand_leaf,
        :nitrogen_deficit,
        :water_deficit,
        # Written at organ temperature by whichever of the three exposure
        # kernels the run enables, and consumed once per day by the
        # reproductive sink.
        :heat_exposure_hours,
        # The same quantity at the lower grain-filling threshold, consumed by
        # terminal heat. Two fields because the thresholds are different
        # physiology, not a tolerance.
        :filling_exposure_hours,
        # 1 when the harvest index, and not the carbon mass cap or the grain
        # already deposited, set storage carbon that day. Diagnostic only.
        :harvest_index_binding,
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
        # Prognostic, monotone, reset at sowing: heat-driven loss of grain set.
        :grain_set_fraction,
        # The same contract for grain FILLING, which terminal heat reduces.
        # Separate state because the harvest index multiplies the two, so
        # neither can stand in for the other.
        :grain_fill_fraction,
        # The third of the same kind, and the only one that is not a fraction:
        # millimetres of rainfall excess the season has accumulated, converted
        # into a recovery loss at harvest once it passes the crop's tolerance.
        # Separate state because lodging, sprouting and harvest loss take a crop
        # that already set and filled its grain.
        :heavy_rain_excess,
        # The fourth, and a SECOND trigger on the same recovery fraction rather
        # than a variant of the third: season-accumulated wind lodging pressure.
        # Separate state because wind years and rain years are, measurably,
        # different years - the per-cell interannual correlation between
        # high-wind and heavy-rain day counts has a median of +0.033 over 886
        # million cropland cell-days - so `heavy_rain_excess` cannot stand in.
        :lodging_exposure,
        # NPP accumulated inside the critical window, which a saturating response
        # turns into a grain number - the one field here that is written
        # over a window and read for the rest of the season rather than reset or
        # monotonically damaged. It replaces a prescribed harvest index that
        # measurement showed to be a constant (`docs/34`), and it is what makes
        # `grain_set_fraction` and `grain_fill_fraction` act on a NUMBER and a
        # WEIGHT respectively instead of both scaling the same constant.
        :window_assimilate,
    )
    @test propertynames(crop.auxiliary.phenology) == (:phu, :winter_type, :fphu)
    @test propertynames(crop.auxiliary.calendar) == (:sowing_date, :prescribed_sowing_date)
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
    @test restart.process_memory.calendar.prescribed_sowing_date ===
          crop.auxiliary.calendar.prescribed_sowing_date
    @test :workspace ∉ propertynames(restart)
    @test :workspace ∉ fieldnames(typeof(output))
    @test propertynames(output.annual) == (
        :yield, :harvest_date,
        :season_gpp, :season_lai_days, :season_length,
        :season_water_deficit, :season_evapotranspiration,
        :harvest_aboveground_carbon,
        # The two diagnostics docs/16 turns on: which side of the flowering
        # window the assimilate came from, and whether the harvest index was the
        # binding constraint at all.
        :window_npp, :hi_binding_days,
        :active_gpp, :active_lai_days, :active_length,
        :active_water_deficit, :active_evapotranspiration,
        :active_window_npp, :active_hi_binding_days,
    )
    @test isempty(fieldnames(typeof(crop.workspace)))
end
