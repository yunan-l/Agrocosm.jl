# C3 keeps the existing parameter-AD solve. C4 uses its own assimilation law;
# the production photosynthesis and canopy kernels are dispatched separately.
@inline _enzyme_weather_lambda(::Val{:C3}, args...) = _enzyme_smooth_lambda_c3(args...)

@inline function _enzyme_weather_lambda(
    ::Val{:C4}, fac::T, vcmax, stress, b, co2, temperature, apar, daylength,
    lpjmlparams, photoparams, upper_bound, iterations, constrain,
) where {T}
    lambda = Agrocosm.compute_lambda_c4_solution(
        fac, vcmax, stress, b, temperature, apar, daylength,
        lpjmlparams, photoparams, upper_bound, iterations,
    )
    for _ in 1:8
        assimilation = Agrocosm.c4_adtmm_scalar_impl(
            lambda, vcmax, stress, b, temperature, apar, daylength,
            lpjmlparams, photoparams,
        )
        phipi = min(one(T), lambda / T(photoparams.lambdamc4))
        light = stress * T(lpjmlparams.alphac4) * apar *
            T(photoparams.cmass) * T(photoparams.cq) / daylength
        je, jc = phipi * light, vcmax / T(24)
        je_slope = lambda < T(photoparams.lambdamc4) ?
            light / T(photoparams.lambdamc4) : zero(T)
        total = je + jc
        root = sqrt(max(zero(T), total * total - T(4) * T(lpjmlparams.theta) * je * jc))
        root_slope = root > zero(T) ?
            (total - T(2) * T(lpjmlparams.theta) * jc) * je_slope / root : zero(T)
        gross_slope = (je_slope - root_slope) * daylength / (T(2) * T(lpjmlparams.theta))
        scale = (temperature + T(273.15)) / T(photoparams.p) * T(8.314) /
            T(photoparams.cmass) * T(1000)
        slope = -fac - (assimilation > zero(T) ? gross_slope * scale : zero(T))
        updated = lambda - (fac * (one(T) - lambda) - assimilation) / slope
        lambda = constrain ? clamp(updated, zero(T), upper_bound) : updated
    end
    return lambda
end

function _check_weather_case(forcing, state, cft, climate, days, harvest_day)
    ndims(forcing) == 3 && size(forcing, 2) == 1 && size(forcing, 3) == 5 ||
        throw(DimensionMismatch("weather attribution requires (day, 1 cell, 5 variables)"))
    length(Agrocosm.crop_prognostic(state).canopy.lai) == 1 ||
        throw(ArgumentError("weather attribution currently supports one cell"))
    Agrocosm.crop_prognostic(state).phenology.is_growing[1] != 0 ||
        throw(ArgumentError("weather attribution starts from an already-established crop"))
    isempty(days) && throw(ArgumentError("growth days must not be empty"))
    first(days) >= 1 && harvest_day == last(days) + 1 && harvest_day <= size(forcing, 1) ||
        throw(ArgumentError("growth days must end immediately before an in-range harvest day"))
    size(climate.temp) == size(forcing)[1:2] || throw(DimensionMismatch("climate shape mismatch"))
    cft.path in (1, 2) || throw(ArgumentError("only C3 and C4 crops are supported"))
    all(isfinite, forcing) || throw(ArgumentError("weather contains non-finite values"))
    for variable in (2, 3, 5)
        all(>=(zero(eltype(forcing))), view(forcing, :, :, variable)) ||
            throw(ArgumentError("precipitation, shortwave and wind must be nonnegative"))
    end
    return nothing
end

function _controlled_climate(climate, forcing)
    return merge(climate, (
        temp = view(forcing, :, :, 1), prec = view(forcing, :, :, 2),
        sw = view(forcing, :, :, 3), lw = view(forcing, :, :, 4),
        wind = view(forcing, :, :, 5),
    ))
end

"""
Replay an already-established single crop through a specified final harvest
using the ordinary lifecycle, including failure and calendar events. Returned
yield is harvested storage C converted to t dry matter/ha. No input is mutated.
The event schedule must be checked before interpreting a fixed-event gradient.
"""
function Agrocosm.weather_harvest_replay(
    forcing::AbstractArray{T, 3}, initial_state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters, parameters::Agrocosm.ModelParameters, climate,
    days::UnitRange{Int}, harvest_day::Int;
    irrigation::Bool = false, nitrogen_limit_vcmax::Bool = true,
    crop_resp_fix::Bool = true,
    replay_end_day::Int = harvest_day,
    diurnal_config = nothing,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false,
    terminal_heat::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
) where {T <: AbstractFloat}
    _check_weather_case(forcing, initial_state, cft, climate, days, harvest_day)
    # Exactly one writer for the exposure fields, the same rule the production
    # constructor enforces. Without it a replay could disagree with the run it
    # is meant to reproduce.
    count((diurnal_config !== nothing, heat_exposure_config !== nothing,
           daily_statistic_exposure)) <= 1 || throw(ArgumentError(
        "sub-daily photosynthesis, the standalone exposure pass and the " *
        "daily-statistic closed form all write heat_exposure_hours; enable at " *
        "most one",
    ))
    organ_temperature && diurnal_config === nothing && heat_exposure_config === nothing &&
        throw(ArgumentError(
            "organ temperature requires sub-daily photosynthesis or the " *
            "standalone heat-exposure pass to be enabled",
        ))
    if organ_temperature
        for field in (:specific_humidity, :surface_pressure)
            hasproperty(climate, field) || throw(ArgumentError(
                "organ temperature requires a `$field` climate field",
            ))
        end
    end
    (diurnal_config === nothing && heat_exposure_config === nothing &&
     !daily_statistic_exposure && !anthesis_heat && !cold_sterility) ||
        hasproperty(climate, :diurnal_range) ||
        throw(ArgumentError("sub-daily integration requires a `diurnal_range` climate field (tasmax - tasmin)"))
    harvest_day <= replay_end_day <= size(forcing, 1) || error("invalid replay end day")
    state = deepcopy(initial_state)
    driver = cft.path == 1 ? Agrocosm.daily_crop_C3! : Agrocosm.daily_crop_C4!
    processes = Agrocosm.ProcessModules(cft, parameters)
    controlled = _controlled_climate(climate, forcing)
    daily = (day = Int[], harvest_event = Int8[], gpp = T[], npp = T[],
        respiration = T[], biological_fixation_cost = T[], et = T[],
        transpiration = T[], soil_evaporation = T[], lai = T[], apar = T[], fphu = T[],
        biomass_carbon = T[], leaf_carbon = T[], root_carbon = T[], mobile_carbon = T[],
        storage_carbon = T[], soil_water = T[], rootzone_available_water = T[],
        water_sufficiency = T[], canopy_conductance = T[], temperature_stress = T[],
        topsoil_temperature = T[], soil_nitrate = T[], soil_ammonium = T[],
        nitrogen_uptake = T[], fertilizer_input = T[], manure_input = T[],
        mineralization = T[], leaching = T[], volatilization = T[],
        leaf_nitrogen = T[], potential_vcmax = T[], vcmax = T[], nitrogen_limitation = T[])
    harvest_days = Int[]
    crop_yield = T(NaN)
    failed = false
    for day in first(days):replay_end_day
        driver(day, day, processes, controlled, state;
            fertilizer = :yes, manure = true, with_tillage = true,
            irrigation, nitrogen_limit_vcmax, crop_resp_fix, diurnal_config,
            organ_temperature, reproductive_sink,
            heat_exposure_config, daily_statistic_exposure, terminal_heat,
            anthesis_heat, cold_sterility,
            update_vernalization_requirement = false, reuse_output = true)
        crop = Agrocosm.crop_prognostic(state)
        fluxes = Agrocosm.crop_fluxes(state)
        photosynthesis = Agrocosm.crop_photosynthesis_auxiliary(state)
        canopy = Agrocosm.crop_canopy_auxiliary(state)
        soil_nitrogen = Agrocosm.soil_nitrogen_prognostic(state)
        nitrogen_fluxes = Agrocosm.soil_nitrogen_fluxes(state)
        # Read-only CPU single-cell diagnostic boundary, outside the AD and
        # production kernels. Auxiliaries retain their daily process-stage value;
        # prognostic stocks are end-of-day. Harvest reset days must be masked
        # when comparing growing-crop mechanisms across different calendars.
        push!(daily.day, day)
        push!(daily.harvest_event, state.output.calendar.harvest_event[1, 1] != 0)
        push!(daily.gpp, fluxes.carbon.gross_assimilation[1])
        push!(daily.npp, fluxes.carbon.npp[1])
        push!(daily.respiration, fluxes.carbon.respiration[1] + fluxes.carbon.leaf_respiration[1])
        push!(daily.biological_fixation_cost, fluxes.carbon.biological_fixation_cost[1])
        push!(daily.et, _daily_transition_observable(state, :et))
        push!(daily.transpiration, sum(fluxes.water.transpiration_layer))
        push!(daily.soil_evaporation, sum(Agrocosm.soil_water_fluxes(state).evaporation))
        push!(daily.lai, canopy.actual_lai[1])
        push!(daily.apar, canopy.apar[1])
        push!(daily.fphu, Agrocosm.crop_phenology_auxiliary(state).fphu[1])
        push!(daily.biomass_carbon, crop.carbon.biomass[1])
        push!(daily.leaf_carbon, crop.carbon.leaf[1])
        push!(daily.root_carbon, crop.carbon.root[1])
        push!(daily.mobile_carbon, crop.carbon.pool[1])
        push!(daily.storage_carbon, crop.carbon.storage[1])
        push!(daily.soil_water, sum(Agrocosm.soil_water_prognostic(state).storage))
        push!(daily.rootzone_available_water, Agrocosm.crop_root_auxiliary(state).zone_available_water[1])
        push!(daily.water_sufficiency, crop.water.sufficiency[1])
        push!(daily.canopy_conductance, canopy.canopy_conductance[1])
        push!(daily.temperature_stress, photosynthesis.temperature_stress[1])
        push!(daily.topsoil_temperature, Agrocosm.soil_thermal_prognostic(state).temperature[1, 1])
        push!(daily.soil_nitrate, sum(soil_nitrogen.nitrate))
        push!(daily.soil_ammonium, sum(soil_nitrogen.ammonium))
        push!(daily.nitrogen_uptake, fluxes.nitrogen.uptake[1])
        push!(daily.fertilizer_input, fluxes.nitrogen.prescribed_fertilizer_input[1])
        push!(daily.manure_input, fluxes.nitrogen.prescribed_manure_input[1])
        push!(daily.mineralization, sum(nitrogen_fluxes.mineralization))
        push!(daily.leaching, nitrogen_fluxes.leaching[1])
        push!(daily.volatilization, nitrogen_fluxes.volatilization[1])
        push!(daily.leaf_nitrogen, crop.nitrogen.leaf[1])
        push!(daily.potential_vcmax, photosynthesis.potential_vcmax[1])
        push!(daily.vcmax, photosynthesis.vcmax[1])
        push!(daily.nitrogen_limitation, photosynthesis.nitrogen_limitation[1])
        # The later failure-termination kernel overwrites the transient event
        # flag even after a normal harvest. The daily calendar output retains
        # both kinds of event and is the production lifecycle record.
        if state.output.calendar.harvest_event[1, 1] != 0
            push!(harvest_days, day)
            crop_yield = Agrocosm.crop_fluxes(state).carbon.yield[1] / T(0.45) * T(0.01)
            failed = Agrocosm.crop_events(state).harvest[1] != 0
            break
        end
    end
    schedule_matches = harvest_days == [harvest_day]
    return (
        yield = crop_yield, harvest_days, schedule_matches, failed, daily, state,
    )
end

function _weather_yield_block(
    forcing, state, cft, parameters, climate, days, layer_depth, terminal::Bool,
    irrigation::Bool, nitrogen_limit_vcmax::Bool, crop_resp_fix::Bool, pathway,
    diurnal_config = nothing, organ_temperature::Bool = false,
    reproductive_sink::Bool = false, heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false, terminal_heat::Bool = false,
    anthesis_heat::Bool = false, cold_sterility::Bool = false,
)
    T = eltype(forcing)
    # A saved post-sowing state still carries the one-day event. Production's
    # cultivate! clears it on the next day. Fixed-event AD skips cultivate!,
    # so clear this marker explicitly or tillage/litter mixing repeats daily.
    fill!(Agrocosm.crop_events(state).sowing, 0)
    for day in days
        _enzyme_continuous_transition!(
            state, cft, parameters, climate, day, :gpp, layer_depth,
            irrigation, nitrogen_limit_vcmax, crop_resp_fix, nitrogen_limit_vcmax,
            forcing, pathway, diurnal_config, organ_temperature, reproductive_sink,
            heat_exposure_config, daily_statistic_exposure, terminal_heat,
            anthesis_heat, cold_sterility,
        )
    end
    # Production harvest_state_kernel! transfers storage carbon directly to
    # harvested yield before any harvest-day allocation. This terminal seed is
    # valid only for the independently verified, fixed calendar harvest.
    return terminal ? Agrocosm.crop_prognostic(state).carbon.storage[1] / T(0.45) * T(0.01) : zero(T)
end

function _weather_reference(forcing, state, cft, parameters, climate, days, harvest_day, irrigation, nitrogen, respiration, diurnal_config = nothing, organ_temperature::Bool = false, reproductive_sink::Bool = false, heat_exposure_config = nothing, daily_statistic_exposure::Bool = false, terminal_heat::Bool = false, anthesis_heat::Bool = false, cold_sterility::Bool = false)
    reference = Agrocosm.weather_harvest_replay(
        forcing, state, cft, parameters, climate, days, harvest_day;
        irrigation, nitrogen_limit_vcmax = nitrogen, crop_resp_fix = respiration, diurnal_config,
        organ_temperature, reproductive_sink,
        heat_exposure_config, daily_statistic_exposure, terminal_heat,
        anthesis_heat, cold_sterility,
    )
    reference.schedule_matches || throw(ArgumentError(
        "fixed-event attribution requires exactly one harvest on day $harvest_day; observed $(reference.harvest_days)",
    ))
    reference.failed && throw(ArgumentError("crop failure is not a differentiable calendar harvest"))
    return reference
end

"""
Compute a blockwise reverse weather gradient of one fixed-calendar harvest.
Parameters, management, CO₂, deposition, and pre-event state are held fixed.
All ordinary-production harvest events are checked first. Event-date or crop
failure discontinuities require separate ordinary counterfactual simulations.
"""
function Agrocosm.enzyme_weather_harvest_gradient(
    forcing::Array{T, 3}, initial_state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters, parameters::Agrocosm.ModelParameters, climate,
    days::UnitRange{Int}, harvest_day::Int;
    block_days::Int = 30, irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = true, crop_resp_fix::Bool = true,
    primal_rtol::Real = 1e-3, primal_atol::Real = 1e-5,
    diurnal_config = nothing,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false,
    terminal_heat::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
) where {T <: AbstractFloat}
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    reference = _weather_reference(forcing, initial_state, cft, parameters, climate,
        days, harvest_day, irrigation, nitrogen_limit_vcmax, crop_resp_fix, diurnal_config,
        organ_temperature, reproductive_sink, heat_exposure_config,
        daily_statistic_exposure, terminal_heat, anthesis_heat, cold_sterility)
    state = deepcopy(initial_state)
    Agrocosm.enzyme_prepare_daily_state!(state)
    layer_depth = Tuple(state.inputs.soil.properties.layer_depth)
    pathway = cft.path == 1 ? Val(:C3) : Val(:C4)
    ranges = [day:min(day + block_days - 1, last(days)) for day in first(days):block_days:last(days)]
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    primal = zero(T)
    for index in eachindex(ranges)
        snapshots[index] = deepcopy(state)
        primal += _weather_yield_block(forcing, state, cft, parameters, climate,
            ranges[index], layer_depth, index == length(ranges), irrigation,
            nitrogen_limit_vcmax, crop_resp_fix, pathway, diurnal_config,
            organ_temperature, reproductive_sink, heat_exposure_config,
            daily_statistic_exposure, terminal_heat, anthesis_heat, cold_sterility)
    end
    isapprox(primal, reference.yield; rtol = primal_rtol, atol = primal_atol) ||
        throw(ArgumentError("weather AD primal $primal differs from production harvest $(reference.yield)"))
    gradient = zeros(T, size(forcing))
    shadow = Agrocosm.enzyme_zero_tangent(state)
    reverse_primal = zero(T)
    for index in reverse(eachindex(ranges))
        block_state = deepcopy(snapshots[index])
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal), _weather_yield_block,
            Enzyme.Duplicated(forcing, gradient), Enzyme.Duplicated(block_state, shadow),
            Enzyme.Const(cft), Enzyme.Const(parameters), Enzyme.Const(climate),
            Enzyme.Const(ranges[index]), Enzyme.Const(layer_depth),
            Enzyme.Const(index == length(ranges)), Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax), Enzyme.Const(crop_resp_fix), Enzyme.Const(pathway),
            Enzyme.Const(diurnal_config), Enzyme.Const(organ_temperature),
            Enzyme.Const(reproductive_sink), Enzyme.Const(heat_exposure_config),
            Enzyme.Const(daily_statistic_exposure), Enzyme.Const(terminal_heat),
            Enzyme.Const(anthesis_heat), Enzyme.Const(cold_sterility),
        )
        reverse_primal += result[2]
    end
    all(isfinite, gradient) || error("weather gradient contains non-finite values")
    return (; primal, reverse_primal, production_yield = reference.yield, gradient,
        variables = Agrocosm.WEATHER_VARIABLES, block_ranges = ranges,
        harvest_day, attribution_scope = :fixed_harvest_conditional)
end

"""One yield block with the CFT rebuilt from `theta` at every block.

The rebuild has to happen INSIDE the differentiated function, not before it, or
`theta` is not on the tape and its gradient is structurally zero. Mirrors
`_enzyme_seasonal_loss_block`.
"""
function _parameter_yield_block(
    theta, parameter_names, forcing, state, base_cft, parameters, climate, days,
    layer_depth, terminal::Bool, irrigation::Bool, nitrogen_limit_vcmax::Bool,
    crop_resp_fix::Bool, pathway, diurnal_config, organ_temperature::Bool,
    reproductive_sink::Bool, heat_exposure_config, daily_statistic_exposure::Bool,
    terminal_heat::Bool, anthesis_heat::Bool, cold_sterility::Bool,
)
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    return _weather_yield_block(
        forcing, state, cft, parameters, climate, days, layer_depth, terminal,
        irrigation, nitrogen_limit_vcmax, crop_resp_fix, pathway, diurnal_config,
        organ_temperature, reproductive_sink, heat_exposure_config,
        daily_statistic_exposure, terminal_heat, anthesis_heat, cold_sterility,
    )
end

"""
Reverse gradient of one fixed-calendar harvest with respect to PROCESS
PARAMETERS rather than weather.

This is the quantity the paper's process-attribution figure needs:
`d(yield)/d(theta)` for the parameters that define each stress mechanism, on the
same fixed-harvest contract as the weather gradient. Weather, management, CO2
and pre-event state are held fixed; `theta` supplies the named `CFTParameters`
fields, rebuilt inside the differentiated block so the tape carries them.

`parameter_names` is a `Tuple{Symbol}` of `CFTParameters` fields. The sink and
terminal-heat parameters live in the tail of a 672-byte struct, past Enzyme's
default type-analysis offset limit, where a reverse derivative comes back as an
exact zero with no error - `AgrocosmEnzymeExt.__init__` raises the limit and
`test/ad/test_cft_offset_gradient.jl` is the guard. A gradient of exactly zero
here should be treated as suspect until that test is known to pass.

Blockwise for the same reason as the weather gradient: the tape for a full
season does not fit. Block gradients ACCUMULATE into one `theta`-shaped vector,
which is correct because every block sees the same parameters - unlike the
weather gradient, where each block owns a disjoint slice of the forcing.
"""
function Agrocosm.enzyme_process_parameter_gradient(
    theta::Vector{T}, parameter_names::Tuple,
    forcing::Array{T, 3}, initial_state::Agrocosm.ModelState,
    base_cft::Agrocosm.CFTParameters, parameters::Agrocosm.ModelParameters, climate,
    days::UnitRange{Int}, harvest_day::Int;
    block_days::Int = 30, irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = true, crop_resp_fix::Bool = true,
    primal_rtol::Real = 1e-3, primal_atol::Real = 1e-5,
    diurnal_config = nothing,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false,
    terminal_heat::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
) where {T <: AbstractFloat}
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta has $(length(theta)) entries for $(length(parameter_names)) names",
    ))
    fields = fieldnames(Agrocosm.CFTParameters)
    for name in parameter_names
        name in fields || throw(ArgumentError("$name is not a CFTParameters field"))
    end
    allunique(parameter_names) || throw(ArgumentError(
        "parameter_names must be unique; a repeated name would sum two " *
        "gradients into one slot",
    ))
    # The reference runs with theta already substituted, so the primal check
    # compares like with like.
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    reference = _weather_reference(forcing, initial_state, cft, parameters, climate,
        days, harvest_day, irrigation, nitrogen_limit_vcmax, crop_resp_fix,
        diurnal_config, organ_temperature, reproductive_sink, heat_exposure_config,
        daily_statistic_exposure, terminal_heat, anthesis_heat, cold_sterility)
    state = deepcopy(initial_state)
    Agrocosm.enzyme_prepare_daily_state!(state)
    layer_depth = Tuple(state.inputs.soil.properties.layer_depth)
    pathway = base_cft.path == 1 ? Val(:C3) : Val(:C4)
    ranges = [day:min(day + block_days - 1, last(days)) for day in first(days):block_days:last(days)]
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    primal = zero(T)
    for index in eachindex(ranges)
        snapshots[index] = deepcopy(state)
        primal += _parameter_yield_block(theta, parameter_names, forcing, state,
            base_cft, parameters, climate, ranges[index], layer_depth,
            index == length(ranges), irrigation, nitrogen_limit_vcmax,
            crop_resp_fix, pathway, diurnal_config, organ_temperature,
            reproductive_sink, heat_exposure_config, daily_statistic_exposure,
            terminal_heat, anthesis_heat, cold_sterility)
    end
    isapprox(primal, reference.yield; rtol = primal_rtol, atol = primal_atol) ||
        throw(ArgumentError("parameter AD primal $primal differs from production harvest $(reference.yield)"))
    gradient = zeros(T, length(theta))
    shadow = Agrocosm.enzyme_zero_tangent(state)
    reverse_primal = zero(T)
    for index in reverse(eachindex(ranges))
        block_state = deepcopy(snapshots[index])
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal), _parameter_yield_block,
            Enzyme.Duplicated(theta, gradient), Enzyme.Const(parameter_names),
            Enzyme.Const(forcing), Enzyme.Duplicated(block_state, shadow),
            Enzyme.Const(base_cft), Enzyme.Const(parameters), Enzyme.Const(climate),
            Enzyme.Const(ranges[index]), Enzyme.Const(layer_depth),
            Enzyme.Const(index == length(ranges)), Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax), Enzyme.Const(crop_resp_fix),
            Enzyme.Const(pathway), Enzyme.Const(diurnal_config),
            Enzyme.Const(organ_temperature), Enzyme.Const(reproductive_sink),
            Enzyme.Const(heat_exposure_config),
            Enzyme.Const(daily_statistic_exposure), Enzyme.Const(terminal_heat),
            Enzyme.Const(anthesis_heat), Enzyme.Const(cold_sterility),
        )
        reverse_primal += result[2]
    end
    all(isfinite, gradient) || error("parameter gradient contains non-finite values")
    return (; primal, reverse_primal, production_yield = reference.yield, gradient,
        parameter_names, theta = copy(theta), block_ranges = ranges, harvest_day,
        attribution_scope = :fixed_harvest_conditional)
end

function Agrocosm.enzyme_weather_forward_directional(
    forcing::Array{T, 3}, direction::Array{T, 3}, initial_state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters, parameters::Agrocosm.ModelParameters, climate,
    days::UnitRange{Int}, harvest_day::Int;
    irrigation::Bool = false, nitrogen_limit_vcmax::Bool = true, crop_resp_fix::Bool = true,
    diurnal_config = nothing,
    organ_temperature::Bool = false,
    reproductive_sink::Bool = false,
    heat_exposure_config = nothing,
    daily_statistic_exposure::Bool = false,
    terminal_heat::Bool = false,
    anthesis_heat::Bool = false,
    cold_sterility::Bool = false,
) where {T <: AbstractFloat}
    size(forcing) == size(direction) || throw(DimensionMismatch("direction shape mismatch"))
    _weather_reference(forcing, initial_state, cft, parameters, climate,
        days, harvest_day, irrigation, nitrogen_limit_vcmax, crop_resp_fix, diurnal_config,
        organ_temperature, reproductive_sink, heat_exposure_config,
        daily_statistic_exposure, terminal_heat, anthesis_heat, cold_sterility)
    state = deepcopy(initial_state)
    Agrocosm.enzyme_prepare_daily_state!(state)
    shadow = Agrocosm.enzyme_zero_tangent(state)
    result = Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ForwardWithPrimal), _weather_yield_block,
        Enzyme.Duplicated, Enzyme.Duplicated(forcing, direction),
        Enzyme.Duplicated(state, shadow), Enzyme.Const(cft), Enzyme.Const(parameters),
        Enzyme.Const(climate), Enzyme.Const(days),
        Enzyme.Const(Tuple(state.inputs.soil.properties.layer_depth)), Enzyme.Const(true),
        Enzyme.Const(irrigation), Enzyme.Const(nitrogen_limit_vcmax),
        Enzyme.Const(crop_resp_fix), Enzyme.Const(cft.path == 1 ? Val(:C3) : Val(:C4)),
        Enzyme.Const(diurnal_config), Enzyme.Const(organ_temperature),
        Enzyme.Const(reproductive_sink), Enzyme.Const(heat_exposure_config),
        Enzyme.Const(daily_statistic_exposure), Enzyme.Const(terminal_heat),
        Enzyme.Const(anthesis_heat), Enzyme.Const(cold_sterility),
    )
    return (primal = result[2], directional = result[1])
end
