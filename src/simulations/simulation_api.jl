"""
    CropSimulation

High-level simulation container. It groups model state, diagnostics, numerical
precision, backend, and run options while preserving direct access to state
objects such as `simulation.crop`, `simulation.soil`, and `simulation.output`.
"""
mutable struct CropSimulation{P, S, D, M, C}
    pft::P
    state::S
    diagnostics::D
    model_parameters::M
    config::C
    simulated_days::Int
end

const _SIMULATION_STATE_PROPERTIES = (
    :climbuf, :crop, :pet, :soil, :managed_land, :daily_weather, :output,
)
const _SIMULATION_DIAGNOSTIC_PROPERTIES = (
    :water_balance, :nitrogen_balance, :carbon_balance, :thermal_balance,
)

function Base.getproperty(simulation::CropSimulation, name::Symbol)
    if name in _SIMULATION_STATE_PROPERTIES
        return getproperty(getfield(simulation, :state), name)
    elseif name in _SIMULATION_DIAGNOSTIC_PROPERTIES
        return getproperty(getfield(simulation, :diagnostics), name)
    end
    return getfield(simulation, name)
end

function Base.propertynames(::CropSimulation, private::Bool = false)
    public = (
        :pft, :model_parameters, :config, :simulated_days,
        _SIMULATION_STATE_PROPERTIES..., _SIMULATION_DIAGNOSTIC_PROPERTIES...,
    )
    return private ? (public..., :state, :diagnostics) : public
end

function _prepare_initial_data(data::NamedTuple, indices, device, T)
    if hasproperty(data, :ModelState) && hasproperty(data, :soilparams)
        return data
    end
    indices === nothing && throw(ArgumentError(
        "indices are required when initialize_simulation receives raw initial data",
    ))
    return InitialDataLoader(data, collect(Int, indices), device; T = T)
end

"""
    initialize_simulation(pft, initial_data; kwargs...)

Create a precision- and backend-consistent crop simulation. `initial_data` may
be either raw input data accepted by `InitialDataLoader` or its normalized
output. Set `diagnostics=false` to avoid allocating daily balance ledgers.
"""
function initialize_simulation(
    pft::PftParameters,
    initial_data::NamedTuple;
    indices = nothing,
    device = identity,
    T::Type{<:AbstractFloat} = Float32,
    days::Integer,
    diagnostics::Bool = true,
    irrigation::Bool = false,
    manure::Bool = false,
    auto_fertilizer::Bool = true,
    nitrogen_limit_vcmax::Bool = false,
    mineral_nitrogen_initialization::Symbol = :lpjml_initsoil,
    c_shift_initialization::Symbol = :lpjml_initsoil,
)
    days > 0 || throw(ArgumentError("days must be positive"))
    prepared = _prepare_initial_data(initial_data, indices, device, T)
    cells = length(prepared.latitude)
    states = init_states!(
        pft, prepared, cells, device;
        T = T,
        mineral_nitrogen_initialization = mineral_nitrogen_initialization,
        c_shift_initialization = c_shift_initialization,
    )
    state = (
        climbuf = states[1],
        crop = states[2],
        pet = states[3],
        soil = states[4],
        managed_land = states[5],
        daily_weather = states[6],
        output = states[7],
    )
    balances = diagnostics ? (
        water_balance = init_water_balance(days, cells, device; T = T),
        nitrogen_balance = init_nitrogen_balance(days, cells, device; T = T),
        carbon_balance = init_carbon_balance(days, cells, device; T = T),
        thermal_balance = init_thermal_balance(days, cells, device; T = T),
    ) : (
        water_balance = nothing,
        nitrogen_balance = nothing,
        carbon_balance = nothing,
        thermal_balance = nothing,
    )
    config = (
        indices = indices === nothing ? nothing : collect(Int, indices),
        device = device,
        T = T,
        days = Int(days),
        irrigation = irrigation,
        manure = manure,
        auto_fertilizer = auto_fertilizer,
        nitrogen_limit_vcmax = nitrogen_limit_vcmax,
    )
    return CropSimulation(
        convert_precision(T, pft), state, balances, ModelParameters(T), config, 0,
    )
end

function _prepare_climate(simulation::CropSimulation, climate::NamedTuple)
    if hasproperty(climate, :swdown) && hasproperty(climate, :lwnet)
        indices = simulation.config.indices
        indices === nothing && throw(ArgumentError(
            "raw climate data require indices in initialize_simulation",
        ))
        return ClimateDataLoader(
            climate, indices, simulation.config.device; T = simulation.config.T,
        )
    end
    return climate
end

"""
    run_simulation!(simulation, climate; start_day=1, end_day=nothing, spinup=true)

Append one continuous climate block. File-local climate rows start at one for
every block, while phenology, outputs, and diagnostics continue from
`simulation.simulated_days`.
"""
function run_simulation!(
    simulation::CropSimulation,
    climate::NamedTuple;
    start_day::Integer = 1,
    end_day::Union{Nothing, Integer} = nothing,
    spinup::Bool = true,
    spinup_years::Integer = 1,
)
    prepared_climate = _prepare_climate(simulation, climate)
    climate_days = size(prepared_climate.temp, 1)
    remaining_days = simulation.config.days - simulation.simulated_days
    remaining_days > 0 || throw(ArgumentError(
        "simulation already contains all $(simulation.config.days) configured days",
    ))
    local_end_day = end_day === nothing ?
        min(climate_days, start_day + remaining_days - 1) : Int(end_day)
    1 <= start_day <= local_end_day <= climate_days || throw(ArgumentError(
        "require 1 <= start_day <= end_day <= $climate_days climate rows",
    ))
    if ndims(prepared_climate.co2) == 1
        required_co2_years = div(local_end_day - 1, 365) + 1
        length(prepared_climate.co2) >= required_co2_years || throw(DimensionMismatch(
            "annual CO₂ forcing has $(length(prepared_climate.co2)) value(s), " *
            "but climate rows through day $local_end_day require $required_co2_years",
        ))
    elseif ndims(prepared_climate.co2) == 2
        size(prepared_climate.co2, 1) >= local_end_day || throw(DimensionMismatch(
            "daily CO₂ forcing has $(size(prepared_climate.co2, 1)) row(s), " *
            "but end_day is $local_end_day",
        ))
    else
        throw(ArgumentError("climate.co2 must be an annual vector or daily matrix"))
    end
    run_days = local_end_day - start_day + 1
    run_days <= remaining_days || throw(DimensionMismatch(
        "only $remaining_days of $(simulation.config.days) configured days remain, requested $run_days",
    ))

    if spinup && simulation.simulated_days == 0 && hasproperty(prepared_climate, :temp_spinup)
        spin_up_climbuf!(
            simulation.pft,
            prepared_climate.temp_spinup,
            simulation.climbuf;
            year_spinup = spinup_years,
        )
    end

    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        auto_fertilizer = simulation.config.auto_fertilizer,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        water_balance = simulation.water_balance,
        nitrogen_balance = simulation.nitrogen_balance,
        carbon_balance = simulation.carbon_balance,
        thermal_balance = simulation.thermal_balance,
        model_parameters = simulation.model_parameters,
        simulation_day_offset = simulation.simulated_days,
        diagnostic_offset = simulation.simulated_days,
    )
    if simulation.pft.path == 1
        daily_crop_C3!(
            start_day, local_end_day, simulation.pft, prepared_climate,
            simulation.climbuf, simulation.crop, simulation.pet, simulation.soil,
            simulation.managed_land, simulation.daily_weather, simulation.output;
            common...,
        )
    elseif simulation.pft.path == 2
        daily_crop_C4!(
            start_day, local_end_day, simulation.pft, prepared_climate,
            simulation.climbuf, simulation.crop, simulation.pet, simulation.soil,
            simulation.managed_land, simulation.daily_weather, simulation.output;
            common...,
        )
    else
        throw(ArgumentError("unsupported photosynthetic pathway $(simulation.pft.path)"))
    end
    simulation.simulated_days += run_days
    return simulation
end

"""Load the `climate` object from one JLD2 file and append it to a simulation."""
function run_simulation!(
    simulation::CropSimulation,
    climate_file::AbstractString;
    kwargs...,
)
    climate = load(climate_file, "climate")
    climate isa NamedTuple || throw(ArgumentError(
        "JLD2 variable `climate` must be a NamedTuple, got $(typeof(climate))",
    ))
    return run_simulation!(simulation, climate; kwargs...)
end

"""
    run_simulation!(simulation, climate_blocks; spinup=true, spinup_years=1)

Run an ordered collection of climate `NamedTuple`s or JLD2 file paths without
concatenating them in memory. The first block may initialize the climate
buffer; subsequent blocks inherit all crop and soil state.
"""
function run_simulation!(
    simulation::CropSimulation,
    climate_blocks::AbstractVector;
    spinup::Bool = true,
    spinup_years::Integer = 1,
)
    isempty(climate_blocks) && throw(ArgumentError("climate_blocks must not be empty"))
    for (block_index, block) in pairs(climate_blocks)
        run_simulation!(
            simulation, block;
            spinup = spinup && block_index == firstindex(climate_blocks),
            spinup_years = spinup_years,
        )
    end
    return simulation
end

_simulation_host(values) = Array(values)

"""Return a compact, backend-independent summary of a completed simulation."""
function simulation_summary(simulation::CropSimulation)
    simulation.simulated_days > 0 || throw(ArgumentError(
        "run the simulation before requesting its summary",
    ))
    npp = _simulation_host(simulation.output.crop.npp)
    lai = _simulation_host(simulation.output.crop.lai)
    biomass = _simulation_host(simulation.output.crop.biomass)
    water = simulation.water_balance
    nitrogen = simulation.nitrogen_balance
    carbon = simulation.carbon_balance
    thermal = simulation.thermal_balance

    balance_summary = water === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(water.residual)),
        cumulative_residual = sum(_simulation_host(water.residual)),
    )
    nitrogen_summary = nitrogen === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(nitrogen.residual)),
        cumulative_residual = sum(_simulation_host(nitrogen.residual)),
        cumulative_leaching_loss = sum(_simulation_host(nitrogen.leaching_loss)),
    )
    carbon_summary = carbon === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(carbon.residual)),
        cumulative_residual = sum(_simulation_host(carbon.residual)),
        cumulative_npp = sum(_simulation_host(carbon.net_primary_production)),
    )
    thermal_summary = thermal === nothing ? nothing : (
        maximum_absolute_energy_residual = maximum(abs, _simulation_host(thermal.energy_residual)),
        maximum_absolute_percolation_energy_residual = maximum(
            abs, _simulation_host(thermal.percolation_energy_residual),
        ),
    )
    return (
        precision = simulation.config.T,
        cells = size(npp, 2),
        simulated_days = simulation.simulated_days,
        crop = (
            cumulative_npp = sum(npp),
            maximum_lai = maximum(lai),
            maximum_biomass = maximum(biomass),
        ),
        water = balance_summary,
        nitrogen = nitrogen_summary,
        carbon = carbon_summary,
        thermal = thermal_summary,
    )
end
