abstract type AbstractExecutionArchitecture end

"""Host execution using ordinary Julia arrays."""
struct HostArchitecture <: AbstractExecutionArchitecture end

"""Array-backed accelerator execution, parameterized by its transfer function."""
struct AcceleratorArchitecture{F} <: AbstractExecutionArchitecture
    device::F
end

"""Compact active-cell domain and its stable external cell identifiers."""
struct ActiveLandDomain{I <: AbstractVector{Int}, C <: AbstractVector{Int32}}
    indices::I
    cell_ids::C

    function ActiveLandDomain(
        indices::AbstractVector{<:Integer},
        cell_ids::AbstractVector{<:Integer} = indices,
    )
        length(indices) == length(cell_ids) || throw(DimensionMismatch(
            "active indices and cell ids must have equal length",
        ))
        isempty(indices) && throw(ArgumentError("active land domain cannot be empty"))
        all(>(0), indices) || throw(ArgumentError("active indices must be positive"))
        allunique(indices) || throw(ArgumentError("active indices must be unique"))
        allunique(cell_ids) || throw(ArgumentError("active cell ids must be unique"))
        return new{Vector{Int}, Vector{Int32}}(Int.(indices), Int32.(cell_ids))
    end
end

"""Precision, architecture, and active-domain contract for one simulation."""
struct ExecutionContext{T <: AbstractFloat, A <: AbstractExecutionArchitecture, D}
    architecture::A
    domain::D
end

function ExecutionContext(
    ::Type{T}, device, indices::AbstractVector{<:Integer}; cell_ids = indices,
) where {T <: AbstractFloat}
    architecture = device === identity ? HostArchitecture() : AcceleratorArchitecture(device)
    domain = ActiveLandDomain(indices, cell_ids)
    return ExecutionContext{T, typeof(architecture), typeof(domain)}(architecture, domain)
end

"""
    SimulationConfiguration

Immutable assembly contract for one `CropSimulation`. It contains only
execution and run choices; scientific parameters remain in `ProcessModules`
and all numerical arrays remain in `ModelState`.
"""
struct SimulationConfiguration{T <: AbstractFloat, D, E}
    indices::Union{Nothing, Vector{Int}}
    device::D
    T::Type{T}
    days::Int
    irrigation::Bool
    manure::Bool
    fertilizer::Symbol
    with_tillage::Bool
    crop_resp_fix::Bool
    nitrogen_limit_vcmax::Bool
    subdaily_photosynthesis::Bool
    subdaily_steps::Int
    diurnal_shape::Symbol
    subdaily_capacity_optimum::Bool
    subdaily_heat_exposure::Bool
    daily_statistic_exposure::Bool
    organ_temperature::Bool
    reproductive_sink::Bool
    anthesis_heat::Bool
    cold_sterility::Bool
    excess_water::Bool
    diurnal_temperature_stress::Bool
    terminal_heat::Bool
    water_sterility::Bool
    water_filling::Bool
    freeze_vernalization_requirement::Bool
    sowing_mode::Symbol
    execution::E
end

function SimulationConfiguration(
    ::Type{T}, device, days::Integer,
    active_indices::AbstractVector{<:Integer}, cell_ids::AbstractVector{<:Integer};
    indices = nothing,
    irrigation::Bool = false,
    manure::Bool = false,
    fertilizer::Symbol = :auto,
    with_tillage::Bool = true,
    crop_resp_fix::Bool = true,
    nitrogen_limit_vcmax::Bool = false,
    subdaily_photosynthesis::Bool = false,
    subdaily_steps::Integer = 24,
    diurnal_shape::Symbol = :sinusoid,
    subdaily_capacity_optimum::Bool = false,
    subdaily_heat_exposure::Bool = false,
    daily_statistic_exposure::Bool = false,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
    excess_water::Bool = false,
    diurnal_temperature_stress::Bool = false,
    terminal_heat::Bool = false,
    water_sterility::Bool = false,
    water_filling::Bool = false,
    freeze_vernalization_requirement::Bool = false,
    sowing_mode::Symbol = :prescribed_sdate,
) where {T <: AbstractFloat}
    days > 0 || throw(ArgumentError("days must be positive"))
    sowing_mode in (:prescribed_sdate, :dynamic_sdate) || throw(ArgumentError(
        "sowing_mode must be :prescribed_sdate or :dynamic_sdate",
    ))
    subdaily_steps >= 1 || throw(ArgumentError("subdaily_steps must be at least 1"))
    diurnal_shape in (:flat, :sinusoid, :daytime_neutral) || throw(ArgumentError(
        "diurnal_shape must be :flat, :sinusoid or :daytime_neutral",
    ))
    # Re-solving Rubisco capacity against the sub-daily light course is a
    # property of that light course, so it is meaningless without it.
    !subdaily_capacity_optimum || subdaily_photosynthesis || throw(ArgumentError(
        "subdaily_capacity_optimum requires subdaily_photosynthesis",
    ))
    # `heat_exposure_hours` must have exactly one writer. The sub-daily
    # assimilation kernels fill it inside the loop they already run; the
    # standalone pass fills it without one. Allowing both would leave the field
    # carrying whichever kernel ran last, which would silently invalidate every
    # ablation cell that reads it - so the combination is rejected rather than
    # ordered.
    writers = count((subdaily_photosynthesis, subdaily_heat_exposure,
                     daily_statistic_exposure))
    writers <= 1 || throw(ArgumentError(
        "subdaily_photosynthesis, subdaily_heat_exposure and " *
        "daily_statistic_exposure all write heat_exposure_hours; enable at " *
        "most one",
    ))
    # The closed form reconstructs the temperature course from the daily mean
    # and range and takes the duration analytically, so there are no sub-steps
    # for a canopy energy balance to be solved at. That is not a limitation to
    # be patched: the cell exists to represent what a model holding only daily
    # aggregates can compute, and giving it leaf temperature would defeat its
    # purpose.
    !(daily_statistic_exposure && organ_temperature) || throw(ArgumentError(
        "daily_statistic_exposure is an air-temperature closed form and has no " *
        "sub-steps to solve a canopy energy balance at; organ_temperature " *
        "requires one of the sub-daily loops",
    ))
    # Leaf temperature is solved per sub-step, so it cannot be switched on by
    # itself - but either sub-daily loop can host it. Inside the assimilation
    # loop it drives the enzyme kinetics as well as the exposure integral;
    # inside the standalone pass it drives the exposure integral alone, which is
    # the configuration that keeps a calibrated daily assimilation kernel while
    # still feeding the sink leaf temperature. Rejecting the combination here
    # keeps the invalid configuration from reaching the kernel, where it could
    # only be ignored.
    !organ_temperature || subdaily_photosynthesis || subdaily_heat_exposure ||
        throw(ArgumentError(
            "organ_temperature requires subdaily_photosynthesis or " *
            "subdaily_heat_exposure",
        ))

    # Sterility is accumulated per sub-step, so the sink needs a sub-daily loop
    # - either the assimilation one or the standalone exposure pass - but not
    # necessarily organ temperature. With organ temperature on it
    # integrates duration at leaf temperature, which is the default and the
    # physically right choice; with it off the same accumulator integrates
    # duration at sub-daily AIR temperature. That second combination is not a
    # mistake to be rejected, it is the ablation cell that separates the sink
    # mechanism from the leaf-air departure that triggers it, and the one an
    # air-temperature-driven GGCM sterility function corresponds to.
    !reproductive_sink || writers >= 1 || throw(ArgumentError(
        "reproductive_sink requires one of subdaily_photosynthesis, " *
        "subdaily_heat_exposure or daily_statistic_exposure to fill " *
        "heat_exposure_hours",
    ))
    # Terminal heat reads `filling_exposure_hours`, which the same three kernels
    # fill in the same loop against a lower threshold, so it has the same
    # prerequisite and no additional one. It is independent of the sink: a run
    # may carry either, both or neither, because grain set and grain filling are
    # separate damage paths and the ablation has to be able to separate them.
    # `water_sterility` and `water_filling` deliberately have NO prerequisite. It reads
    # `water.sufficiency`, which `transpiration!` writes on every day of every
    # configuration, so there is no exposure source to enable first - unlike
    # every other mechanism this project added. Nothing to validate here is the
    # correct outcome, not an omission.
    !terminal_heat || writers >= 1 || throw(ArgumentError(
        "terminal_heat requires one of subdaily_photosynthesis, " *
        "subdaily_heat_exposure or daily_statistic_exposure to fill " *
        "filling_exposure_hours",
    ))
    execution = ExecutionContext(T, device, active_indices; cell_ids)
    source_indices = indices === nothing ? nothing : Int.(indices)
    return SimulationConfiguration{
        T, typeof(device), typeof(execution),
    }(
        source_indices, device, T, Int(days), irrigation, manure, fertilizer,
        with_tillage, crop_resp_fix, nitrogen_limit_vcmax,
        subdaily_photosynthesis, Int(subdaily_steps), diurnal_shape,
        subdaily_capacity_optimum, subdaily_heat_exposure,
        daily_statistic_exposure, organ_temperature, reproductive_sink,
        anthesis_heat, cold_sterility, excess_water, diurnal_temperature_stress,
        terminal_heat, water_sterility, water_filling,
        freeze_vernalization_requirement,
        sowing_mode, execution,
    )
end

"""
    diurnal_configuration(config)

Zero-size `DiurnalConfig` for this run, or `nothing` when sub-daily integration
is switched off. `nothing` routes every assimilation call to the existing daily
kernel, so production is bitwise unchanged.
"""
diurnal_configuration(config::SimulationConfiguration) =
    config.subdaily_photosynthesis ?
        DiurnalConfig(; steps = config.subdaily_steps,
                        shape = diurnal_shape_code(config.diurnal_shape),
                        capacity_optimum = config.subdaily_capacity_optimum) :
        nothing

"""
    daily_statistic_exposure_enabled(config)

Whether the closed-form daily-statistic exposure path fills
`heat_exposure_hours`. A plain `Bool` rather than a `DiurnalConfig`, because the
closed form has no sub-step count and no shape to carry: it is the analytic
duration for the sinusoid `diurnal_temperature` reconstructs.
"""
daily_statistic_exposure_enabled(config::SimulationConfiguration) =
    config.daily_statistic_exposure

"""
    heat_exposure_configuration(config)

Zero-size `DiurnalConfig` for the standalone exposure pass, or `nothing` when it
is switched off. `nothing` makes `heat_exposure!` a no-op, so the driver can
call it unconditionally.

Shares `subdaily_steps` and `diurnal_shape` with the assimilation loop because
they are properties of the run's sub-daily resolution, not of either loop. The
capacity solve is not: it re-solves Rubisco capacity against the light course,
which this pass does not integrate, so it is off here regardless of the run's
setting - and `subdaily_capacity_optimum` requires `subdaily_photosynthesis`,
which is mutually exclusive with this switch, so it is already false in any
configuration that reaches this line.
"""
heat_exposure_configuration(config::SimulationConfiguration) =
    config.subdaily_heat_exposure ?
        DiurnalConfig(; steps = config.subdaily_steps,
                        shape = diurnal_shape_code(config.diurnal_shape),
                        capacity_optimum = false) :
        nothing

float_type(::ExecutionContext{T}) where {T} = T
array_device(::HostArchitecture) = identity
array_device(architecture::AcceleratorArchitecture) = architecture.device
array_device(context::ExecutionContext) = array_device(context.architecture)
architecture_name(::HostArchitecture) = :cpu
architecture_name(::AcceleratorArchitecture) = :accelerator
architecture_name(context::ExecutionContext) = architecture_name(context.architecture)

function Base.show(io::IO, domain::ActiveLandDomain)
    print(io, "ActiveLandDomain(", length(domain.indices), " cells)")
end

function Base.show(io::IO, context::ExecutionContext)
    print(
        io,
        "ExecutionContext(", architecture_name(context), ", ",
        float_type(context), ", ", length(context.domain.indices), " cells)",
    )
end

"""Machine-readable metadata for a model or output variable."""
struct VariableSpec
    path::Tuple{Vararg{Symbol}}
    role::Symbol
    dimensions::Tuple{Vararg{Symbol}}
    units::String
    description::String
end

const _OUTPUT_VARIABLE_METADATA = Dict{Tuple{Symbol, Symbol}, NamedTuple}(
    (:crop, :gpp) => (units = "gC m-2 day-1", description = "Gross primary production"),
    (:crop, :npp) => (units = "gC m-2 day-1", description = "Net primary production"),
    (:crop, :lambda) => (units = "1", description = "Intercellular-to-ambient CO2 ratio"),
    (:crop, :potential_vcmax) => (units = "gC m-2 day-1", description = "Potential maximum carboxylation capacity"),
    (:crop, :vcmax) => (units = "gC m-2 day-1", description = "Realized maximum carboxylation capacity"),
    (:crop, :nitrogen_limitation) => (units = "1", description = "Nitrogen limitation factor"),
    (:crop, :respiration) => (units = "gC m-2 day-1", description = "Plant respiration"),
    (:crop, :biomass) => (units = "gC m-2", description = "Live crop carbon"),
    (:crop, :lai) => (units = "m2 m-2", description = "Leaf area index"),
    (:crop, :storage_carbon) => (units = "gC m-2", description = "Storage-organ carbon"),
    (:crop, :yield) => (units = "gC m-2 year-1", description = "Harvested storage-organ carbon"),
    (:crop, :season_gpp) => (units = "gC m-2", description = "Harvest-season cumulative gross primary production"),
    (:crop, :season_lai_days) => (units = "m2 m-2 day", description = "Harvest-season cumulative leaf area index"),
    (:crop, :season_length) => (units = "day", description = "Active crop days in the harvested season"),
    # Named for LPJmL's `wdf`, but the quantity is SUFFICIENCY: 100 means
    # supply met demand. `compute_water_sufficiency` returns
    # clamp(100 * sum(min(supply, demand)) / sum(demand), 0, 100), so a
    # larger value is a wetter season, not a drier one. Reading these two
    # the way their names read inverts every drought figure.
    (:crop, :season_water_deficit) => (units = "% day", description = "Harvest-season sum of the daily season-cumulative water sufficiency (100 = demand met)"),
    (:crop, :season_evapotranspiration) => (units = "mm", description = "Harvest-season cumulative evapotranspiration"),
    (:crop, :harvest_aboveground_carbon) => (units = "gC m-2", description = "Live above-ground crop carbon immediately before harvest"),
    # NPP over the flowering window, gated not weighted: grain number responds
    # to assimilate supply over the critical period. See docs/16.
    (:crop, :window_npp) => (units = "gC m-2", description = "Net primary production accumulated on growing days inside the flowering window"),
    # Counts the days a mechanism multiplying the harvest index could have had
    # any effect at all; on the other days the carbon mass cap or the grain
    # already deposited set storage, and the index is not the binding constraint.
    (:crop, :hi_binding_days) => (units = "day", description = "Growing days on which the harvest index, not the carbon mass cap or the already-deposited grain, set storage carbon"),
    (:crop, :fphu) => (units = "1", description = "Fraction of potential heat units"),
    (:crop, :water_deficit) => (units = "%", description = "Season-cumulative crop water sufficiency (100 = demand met)"),
    (:crop, :growing_mask) => (units = "1", description = "Active crop-stand mask"),
    (:soil, :ecosystem_respiration) => (
        units = "gC m-2 day-1",
        description = "Ecosystem respiration: plant plus heterotrophic soil respiration",
    ),
    (:soil, :heterotrophic_respiration) => (
        units = "gC m-2 day-1",
        description = "Heterotrophic litter and soil respiration",
    ),
    (:soil, :evapotranspiration) => (
        units = "mm day-1",
        description = "Total land-surface evapotranspiration",
    ),
    (:calendar, :harvesting_mask) => (units = "1", description = "Harvest condition mask"),
    (:calendar, :harvesting_year) => (units = "year", description = "Simulation harvest year"),
    (:calendar, :harvest_date) => (units = "day_of_year", description = "Harvest day of year"),
    (:calendar, :sowing_event) => (units = "1", description = "Daily sowing event"),
    (:calendar, :harvest_event) => (units = "1", description = "Daily harvest event"),
)

function output_variable_spec(group::Symbol, field::Symbol)
    metadata = get(_OUTPUT_VARIABLE_METADATA, (group, field), nothing)
    isnothing(metadata) && throw(ArgumentError("missing metadata for output variable $group.$field"))
    frequency = field in (
        :yield, :season_gpp, :season_lai_days, :season_length,
        :season_water_deficit, :season_evapotranspiration,
        :harvest_aboveground_carbon, :window_npp, :hi_binding_days,
        :harvest_date, :harvesting_year,
    ) ? :annual : :daily
    return VariableSpec(
        (:output, group, field), :output, (:time, :cell),
        metadata.units, metadata.description,
    ), frequency
end

function _state_array_specs!(specs, value, path::Tuple, role::Symbol, cells::Int)
    if value isa AbstractArray
        dimensions = ntuple(
            index -> index == ndims(value) && size(value, index) == cells ?
                :cell : Symbol(:axis_, index),
            ndims(value),
        )
        push!(specs, VariableSpec(path, role, dimensions, "", ""))
        return specs
    end
    value isa NamedTuple || isstructtype(typeof(value)) || return specs
    for name in propertynames(value)
        _state_array_specs!(specs, getproperty(value, name), (path..., name), role, cells)
    end
    return specs
end

"""Return the numerical state inventory grouped by lifecycle role."""
function state_schema(state::ModelState)
    specs = VariableSpec[]
    cells = length(state.inputs.weather.temp)
    for role in (:prognostic, :fluxes, :auxiliary, :inputs, :events, :workspace)
        _state_array_specs!(specs, getproperty(state, role), (role,), role, cells)
    end
    paths = getfield.(specs, :path)
    allunique(paths) || throw(ArgumentError("model state contains duplicate variable paths"))
    return specs
end

"""Validate that every runtime state array ends in the active-cell dimension."""
function validate_state_schema(state::ModelState, cells::Integer)
    cells > 0 || throw(ArgumentError("cells must be positive"))
    for spec in state_schema(state)
        value = foldl(getproperty, spec.path; init = state)
        :cell in spec.dimensions || continue
        size(value, ndims(value)) == cells || throw(DimensionMismatch(
            "$(join(spec.path, '.')) ends in size $(size(value, ndims(value))); expected $cells active cells",
        ))
    end
    return state
end
