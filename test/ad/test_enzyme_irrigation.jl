function _irrigated_daily_transition_objective(
    theta,
    state,
    base_cft,
    global_parameters,
    climate,
    parameter_names,
    day,
    observable,
    layer_depth,
    irrigation_mode,
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
        irrigation = true,
        irrigation_mode,
    )
end

function _set_frozen_irrigated_state!(state)
    water = Agrocosm.soil_water_prognostic(state)
    field_capacity = Agrocosm.soil_water_auxiliary(state).field_capacity
    layer_depth = Agrocosm.soil_properties(state).layer_depth
    target_storage = field_capacity .* reshape(layer_depth, :, 1)
    water.storage .= target_storage .* 0.5f0
    water.ice_storage .= target_storage .* 0.25f0
    Agrocosm.partition_soil_water_ice!(state)
    return target_storage
end

@testset "Enzyme irrigated frozen-soil reset" begin
    template = _daily_transition_fixture()
    extension = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)
    @test extension !== nothing

    production_state = deepcopy(template.state)
    enzyme_state = deepcopy(template.state)
    target_storage = _set_frozen_irrigated_state!(production_state)
    _set_frozen_irrigated_state!(enzyme_state)

    soil_evapotranspiration!(production_state, production_state; irrigation = true)
    extension._enzyme_soil_evapotranspiration!(enzyme_state; irrigation = true)

    production_water = Agrocosm.soil_water_prognostic(production_state)
    enzyme_water = Agrocosm.soil_water_prognostic(enzyme_state)
    @test enzyme_water.storage == production_water.storage
    @test enzyme_water.ice_storage == production_water.ice_storage
    @test enzyme_water.storage + enzyme_water.ice_storage ≈ target_storage

    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    target_storage = _set_frozen_irrigated_state!(state)
    primal = _irrigated_daily_transition_objective(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        :et,
        template.layer_depth,
        :reset,
    )
    water = Agrocosm.soil_water_prognostic(state)
    @test isfinite(primal)
    @test water.storage + water.ice_storage ≈ target_storage

    reverse_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(reverse_state)
    _set_frozen_irrigated_state!(reverse_state)
    gradient = zeros(Float32, length(theta))
    reverse_result = Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        _irrigated_daily_transition_objective,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Duplicated(reverse_state, enzyme_zero_tangent(reverse_state)),
        Enzyme.Const(cft1),
        Enzyme.Const(template.global_parameters),
        Enzyme.Const(template.climate),
        Enzyme.Const(parameter_names),
        Enzyme.Const(template.day),
        Enzyme.Const(:et),
        Enzyme.Const(template.layer_depth),
        Enzyme.Const(:reset),
    )
    @test isfinite(reverse_result[2])
    @test all(isfinite, gradient)

    forward_gradient = enzyme_forward_gradient(
        _irrigated_daily_transition_objective,
        theta,
        () -> begin
            fresh_state = deepcopy(template.state)
            enzyme_prepare_daily_state!(fresh_state)
            _set_frozen_irrigated_state!(fresh_state)
            return fresh_state, enzyme_zero_tangent(fresh_state)
        end,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        :et,
        template.layer_depth,
        :reset,
    )
    @test gradient ≈ forward_gradient rtol = 5.0f-2 atol = 1.0f-4
end

@testset "Enzyme surface irrigation matches production input" begin
    template = _daily_transition_fixture()
    extension = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)
    state = deepcopy(template.state)
    _set_frozen_irrigated_state!(state)
    pedotransfer!(state; lpjmlparams = template.global_parameters.lpjml)
    precipitation = Float32[1]
    canopy_interception = Float32[0]
    applied = zeros(Float32, 1)
    production_input = similar(precipitation)
    Agrocosm.surface_irrigation!(
        applied, production_input, state, precipitation, canopy_interception, cft1,
    )
    enzyme_input = extension._enzyme_surface_irrigation_input(
        state, cft1, precipitation, canopy_interception,
    )
    @test applied[1] > 0
    @test enzyme_input == production_input

    theta = Float32[cft1.gmin, cft1.b]
    parameter_names = (:gmin, :b)
    forward_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(forward_state)
    _set_frozen_irrigated_state!(forward_state)
    primal = enzyme_daily_transition_objective(
        theta, forward_state, cft1, template.global_parameters, template.climate,
        parameter_names, template.day, :et, template.layer_depth;
        irrigation = true,
    )
    @test isfinite(primal)

    production_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(production_state)
    _set_frozen_irrigated_state!(production_state)
    daily_crop_C3!(
        template.day, template.day,
        ProcessModules(cft1, template.global_parameters),
        template.climate, production_state;
        fertilizer = :yes, manure = true, with_tillage = true,
        irrigation = true, update_vernalization_requirement = false,
        reuse_output = true,
    )
    @test production_state.output.soil.irrigation_applied[1, 1] > 0
    @test primal ≈ extension._daily_transition_observable(production_state, :et) rtol = 5.0f-5 atol = 1.5f-4
    @test Agrocosm.soil_water_prognostic(forward_state).storage ≈
        Agrocosm.soil_water_prognostic(production_state).storage rtol = 5.0f-5 atol = 1.5f-4

    reverse_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(reverse_state)
    _set_frozen_irrigated_state!(reverse_state)
    gradient = zeros(Float32, length(theta))
    result = Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        _irrigated_daily_transition_objective,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Duplicated(reverse_state, enzyme_zero_tangent(reverse_state)),
        Enzyme.Const(cft1), Enzyme.Const(template.global_parameters),
        Enzyme.Const(template.climate), Enzyme.Const(parameter_names),
        Enzyme.Const(template.day), Enzyme.Const(:et),
        Enzyme.Const(template.layer_depth), Enzyme.Const(:surface),
    )
    @test isfinite(result[2])
    @test all(isfinite, gradient)
    @test result[2] ≈ primal rtol = 1.0f-5
end
