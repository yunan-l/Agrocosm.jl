using Agrocosm
using Enzyme
using LinearAlgebra
using Test

function _daily_transition_fixture(climate_days::Int = 16, initial_end_day::Int = 10)
    T = Float32
    cells = 1
    layers = 5
    days = climate_days
    initial_data = (
        latitude = T[45],
        soilparams = (
            ph = T[6.5],
            w_sat = fill(T(0.45), layers, cells),
            sand = reshape(T[0.4], 1, cells),
            clay = reshape(T[0.2], 1, cells),
            tdiff_0 = T[0.7],
            tdiff_15 = T[0.75],
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        ModelState = (
            crop = (
                sdate = Int32[1],
                phu = T[543],
                manure = zeros(T, cells),
                fertilizer = T[24.55],
                residuefrac = T[0.67],
            ),
            u0 = (
                swc = reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], layers, cells),
                litc = reshape(T[0.13, 187.5, 225.36], 3, cells),
                fastc = reshape(T[548.97, 368.27, 313.79, 377.55, 344.65], layers, cells),
                slowc = reshape(T[1218.62, 753.33, 660.10, 792.63, 736.38], layers, cells),
                litn = reshape(T[0.0047, 6.47, 9.47], 3, cells),
                fastn = reshape(T[36.60, 24.55, 20.92, 25.17, 22.98], layers, cells),
                slown = reshape(T[81.24, 50.22, 44.01, 52.84, 49.09], layers, cells),
            ),
        ),
    )
    climbuf, crop, pet, soil, managed_land, weather, output = init_states!(
        cft1, initial_data, cells, identity; T,
    )
    climbuf.atemp .= T(10)
    climbuf.temp .= T(10)
    climbuf.atemp_mean .= T(10)
    climate = (
        temp = fill(T(15), days, cells),
        prec = fill(T(1), days, cells),
        sw = fill(T(180), days, cells),
        lw = fill(T(-40), days, cells),
        wind = fill(T(2), days, cells),
        no3_deposition = fill(T(0.01), days, cells),
        nh4_deposition = fill(T(0.02), days, cells),
        co2 = T[400],
    )
    state = model_state(climbuf, crop, pet, soil, managed_land, weather, output)
    global_parameters = ModelParameters(T)
    processes = ProcessModules(cft1, global_parameters)
    if initial_end_day > 0
        daily_crop_C3!(1, initial_end_day, processes, climate, state;
            fertilizer = :yes,
            manure = true,
            with_tillage = true,
            update_vernalization_requirement = false,
            reuse_output = true,
        )
    end
    enzyme_prepare_daily_state!(state)
    layer_depth = Tuple(state.inputs.soil.properties.layer_depth)
    return (; state, climate, global_parameters, layer_depth, day = initial_end_day + 1)
end

function _daily_transition_value(template, theta, parameter_names, observable)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_daily_transition_objective(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        observable,
        template.layer_depth,
    )
end

function _production_daily_transition_value(state, observable)
    crop_flux = state.fluxes.crop
    if observable === :gpp
        return crop_flux.carbon.gross_assimilation[1]
    elseif observable === :reco
        return crop_flux.carbon.respiration[1] +
            crop_flux.carbon.leaf_respiration[1] +
            state.fluxes.soil.carbon.heterotrophic_respiration[1]
    elseif observable === :et
        crop_water = crop_flux.water
        value = crop_water.interception[1] + state.fluxes.soil.surface_litter.evaporation[1]
        value += crop_water.transpiration_layer[1, 1] + state.fluxes.soil.water.evaporation[1, 1]
        value += crop_water.transpiration_layer[2, 1] + state.fluxes.soil.water.evaporation[2, 1]
        value += crop_water.transpiration_layer[3, 1] + state.fluxes.soil.water.evaporation[3, 1]
        value += crop_water.transpiration_layer[4, 1] + state.fluxes.soil.water.evaporation[4, 1]
        value += crop_water.transpiration_layer[5, 1] + state.fluxes.soil.water.evaporation[5, 1]
        return value
    end
    throw(ArgumentError("unsupported observable $observable"))
end

function _multi_day_transition_objective(
    theta,
    state,
    base_cft,
    global_parameters,
    climate,
    parameter_names,
    days::NTuple{3, Int},
    observable,
    layer_depth,
)
    value = zero(eltype(theta))
    value += enzyme_daily_transition_objective(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, days[1], observable, layer_depth,
    )
    value += enzyme_daily_transition_objective(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, days[2], observable, layer_depth,
    )
    value += enzyme_daily_transition_objective(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, days[3], observable, layer_depth,
    )
    return value
end

function _multi_day_transition_value(template, theta, parameter_names, days, observable)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return _multi_day_transition_objective(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        observable,
        template.layer_depth,
    )
end

function _nitrogen_limited_daily_objective(
    theta,
    state,
    base_cft,
    global_parameters,
    climate,
    parameter_names,
    day,
    observable,
    layer_depth,
)
    return enzyme_daily_transition_objective(
        theta,
        state,
        base_cft,
        global_parameters,
        climate,
        parameter_names,
        day,
        observable,
        layer_depth;
        nitrogen_limit_vcmax = true,
        crop_resp_fix = false,
    )
end

function _nitrogen_limited_daily_value(template, theta, parameter_names, observable)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return _nitrogen_limited_daily_objective(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        observable,
        template.layer_depth,
    )
end

function _nitrogen_limited_transition_fixture()
    template = _daily_transition_fixture()
    soil_nitrogen = Agrocosm.soil_nitrogen_prognostic(template.state)
    soil_nitrogen.nitrate .= 1.0f-4
    soil_nitrogen.ammonium .= 1.0f-4
    plant_nitrogen = Agrocosm.crop_prognostic(template.state).nitrogen
    plant_nitrogen.total .= 2.0f0
    plant_nitrogen.leaf .= 0.8f0
    plant_nitrogen.root .= 0.6f0
    plant_nitrogen.pool .= 0.4f0
    plant_nitrogen.storage .= 0.2f0
    return template
end

@testset "Enzyme one-day ModelState transition" begin
    template = _daily_transition_fixture()
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end

    for observable in (:gpp, :reco, :et)
        production_state = deepcopy(template.state)
        enzyme_prepare_daily_state!(production_state)
        production_processes = ProcessModules(cft1, template.global_parameters)
        daily_crop_C3!(
            template.day, template.day, production_processes,
            template.climate, production_state;
            fertilizer = :yes,
            manure = true,
            with_tillage = true,
            update_vernalization_requirement = false,
            reuse_output = true,
        )
        production_value = _production_daily_transition_value(production_state, observable)
        adapter_value = _daily_transition_value(template, theta, parameter_names, observable)
        # The AD-only lambda solve continuously refines the production bisection
        # root, so parity is checked at the production solver tolerance.
        @test adapter_value ≈ production_value rtol = 5.0f-5 atol = 1.5f-4

        ad = enzyme_forward_directional(
            enzyme_daily_transition_objective,
            theta,
            direction,
            state_factory,
            cft1,
            template.global_parameters,
            template.climate,
            parameter_names,
            template.day,
            observable,
            template.layer_depth;
            return_primal = true,
        )
        @test ad.primal ≈ adapter_value rtol = 1.0f-6 atol = 1.0f-7

        gradient = enzyme_forward_gradient(
            enzyme_daily_transition_objective,
            theta,
            state_factory,
            cft1,
            template.global_parameters,
            template.climate,
            parameter_names,
            template.day,
            observable,
            template.layer_depth,
        )
        projection = dot(gradient, direction)
        projection_scale = sum(abs, gradient .* direction)
        @test abs(ad.directional - projection) <=
            2.0f-5 * projection_scale + 1.0f-6

        finite_difference = similar(theta)
        for index in eachindex(theta)
            step = 1.0f-3 * max(abs(theta[index]), 1.0f0)
            delta = zeros(Float32, length(theta))
            delta[index] = step
            finite_difference[index] = (
                _daily_transition_value(template, theta .+ delta, parameter_names, observable) -
                _daily_transition_value(template, theta .- delta, parameter_names, observable)
            ) / (2.0f0 * step)
        end
        @test gradient ≈ finite_difference rtol = 2.0f-2 atol = 1.0f-5
    end

end

@testset "Enzyme crop respiration mode is an explicit constant" begin
    template = _daily_transition_fixture()
    plant_nitrogen = Agrocosm.crop_prognostic(template.state).nitrogen
    plant_nitrogen.root .= 0.0f0
    plant_nitrogen.storage .= 0.0f0
    plant_nitrogen.pool .= 0.0f0
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]

    function enzyme_reco(; kwargs...)
        state = deepcopy(template.state)
        enzyme_prepare_daily_state!(state)
        return enzyme_daily_transition_objective(
            theta,
            state,
            cft1,
            template.global_parameters,
            template.climate,
            parameter_names,
            template.day,
            :reco,
            template.layer_depth;
            kwargs...,
        )
    end

    default_reco = enzyme_reco()
    fixed_reco = enzyme_reco(; crop_resp_fix = true)
    dynamic_reco = enzyme_reco(; crop_resp_fix = false)
    @test default_reco == fixed_reco
    @test abs(fixed_reco - dynamic_reco) > 1.0f-4

    production_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(production_state)
    daily_crop_C3!(
        template.day,
        template.day,
        ProcessModules(cft1, template.global_parameters),
        template.climate,
        production_state;
        fertilizer = :yes,
        manure = true,
        with_tillage = true,
        crop_resp_fix = false,
        update_vernalization_requirement = false,
        reuse_output = true,
    )
    production_reco = _production_daily_transition_value(production_state, :reco)
    @test dynamic_reco ≈ production_reco rtol = 5.0f-5 atol = 1.5f-4
end

@testset "Enzyme multi-day state propagation" begin
    template = _daily_transition_fixture()
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)
    days = (11, 12, 13)
    observable = :gpp

    production_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(production_state)
    production_processes = ProcessModules(cft1, template.global_parameters)
    production_value = zero(Float32)
    for day in days
        daily_crop_C3!(
            day, day, production_processes,
            template.climate, production_state;
            fertilizer = :yes,
            manure = true,
            with_tillage = true,
            update_vernalization_requirement = false,
            reuse_output = true,
        )
        production_value += _production_daily_transition_value(production_state, observable)
    end
    adapter_value = _multi_day_transition_value(template, theta, parameter_names, days, observable)
    # The AD-only Newton lambda solve is continuous, while production retains
    # LPJmL-compatible bisection. With nitrogen-limited Vcmax, that bounded
    # daily numerical difference feeds back through organ nitrogen and can
    # accumulate across days without affecting AD self-consistency checks.
    @test adapter_value ≈ production_value rtol = 2.0f-2 atol = 5.0f-3

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    ad = enzyme_forward_directional(
        _multi_day_transition_objective,
        theta,
        direction,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        observable,
        template.layer_depth;
        return_primal = true,
    )
    @test ad.primal ≈ adapter_value rtol = 1.0f-6 atol = 1.0f-7

    gradient = enzyme_forward_gradient(
        _multi_day_transition_objective,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        observable,
        template.layer_depth,
    )
    @test ad.directional ≈ dot(gradient, direction) rtol = 1.0f-5 atol = 1.0f-6

    finite_difference = similar(theta)
    for index in eachindex(theta)
        step = 1.0f-3 * max(abs(theta[index]), 1.0f0)
        delta = zeros(Float32, length(theta))
        delta[index] = step
        finite_difference[index] = (
            _multi_day_transition_value(template, theta .+ delta, parameter_names, days, observable) -
            _multi_day_transition_value(template, theta .- delta, parameter_names, days, observable)
        ) / (2.0f0 * step)
    end
    @test gradient ≈ finite_difference rtol = 2.0f-2 atol = 1.0f-5
end

@testset "Enzyme explicit nitrogen-limited transition" begin
    template = _nitrogen_limited_transition_fixture()
    parameter_names = (:knstore, :gmin)
    theta = Float32[cft1.knstore, cft1.gmin]
    observable = :gpp

    production_state = deepcopy(template.state)
    enzyme_state = deepcopy(template.state)
    production_nitrogen = Agrocosm.crop_prognostic(production_state).nitrogen
    enzyme_nitrogen = Agrocosm.crop_prognostic(enzyme_state).nitrogen
    production_nitrogen.pending_manure .= 4.0f0
    enzyme_nitrogen.pending_manure .= 4.0f0
    deferred_fertilizer = only(production_nitrogen.pending_fertilizer)
    enzyme_prepare_daily_state!(production_state)
    enzyme_prepare_daily_state!(enzyme_state)
    daily_crop_C3!(
        template.day, template.day,
        ProcessModules(cft1, template.global_parameters),
        template.climate,
        production_state;
        fertilizer = :yes,
        manure = true,
        with_tillage = true,
        nitrogen_limit_vcmax = true,
        crop_resp_fix = false,
        update_vernalization_requirement = false,
        reuse_output = true,
    )
    enzyme_daily_transition_objective(
        theta,
        enzyme_state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        observable,
        template.layer_depth;
        nitrogen_limit_vcmax = true,
        crop_resp_fix = false,
    )
    @test only(production_nitrogen.pending_fertilizer) == 0.0f0
    @test only(enzyme_nitrogen.pending_fertilizer) == 0.0f0
    @test only(production_nitrogen.pending_manure) == 0.0f0
    @test only(enzyme_nitrogen.pending_manure) == 0.0f0
    @test only(Agrocosm.crop_fluxes(production_state).nitrogen.prescribed_fertilizer_input) ≈
        deferred_fertilizer
    @test only(Agrocosm.crop_fluxes(enzyme_state).nitrogen.prescribed_fertilizer_input) ≈
        deferred_fertilizer
    @test only(Agrocosm.crop_fluxes(production_state).nitrogen.prescribed_manure_input) ≈ 4.0f0
    @test only(Agrocosm.crop_fluxes(enzyme_state).nitrogen.prescribed_manure_input) ≈ 4.0f0

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end

    gradient = enzyme_forward_gradient(
        _nitrogen_limited_daily_objective,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        observable,
        template.layer_depth,
    )
    @test all(isfinite, gradient)
    @test any(!iszero, gradient)

    finite_difference = similar(theta)
    for index in eachindex(theta)
        step = 1.0f-3 * max(abs(theta[index]), 1.0f0)
        delta = zeros(Float32, length(theta))
        delta[index] = step
        finite_difference[index] = (
            _nitrogen_limited_daily_value(
                template, theta .+ delta, parameter_names, observable,
            ) -
            _nitrogen_limited_daily_value(
                template, theta .- delta, parameter_names, observable,
            )
        ) / (2.0f0 * step)
    end
    @test gradient ≈ finite_difference rtol = 3.0f-2 atol = 2.0f-5
end
