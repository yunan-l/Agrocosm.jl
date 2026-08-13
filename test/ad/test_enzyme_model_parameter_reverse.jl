const _AgrocosmEnzymeExtension = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)

function _model_parameter_replacement_sum(theta, base_parameters, parameter_names)
    parameters = _AgrocosmEnzymeExtension._replace_model_parameters(
        base_parameters, theta, parameter_names,
    )
    lpjml = parameters.lpjml
    return lpjml.k_soil10.fast + lpjml.PRIESTLEY_TAYLOR +
        lpjml.soildepth_evap + lpjml.soil_infil
end

@testset "Enzyme model-parameter replacement remains differentiable" begin
    parameter_names = (:k_soil10_fast, :PRIESTLEY_TAYLOR, :soildepth_evap, :soil_infil)
    base_parameters = ModelParameters(Float32)
    lpjml = base_parameters.lpjml
    theta = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    gradient = zeros(Float32, length(theta))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        _model_parameter_replacement_sum,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Const(base_parameters),
        Enzyme.Const(parameter_names),
    )

    @test gradient == ones(Float32, length(theta))
end

@testset "Enzyme blockwise soil reverse" begin
    template = _daily_transition_fixture(14, 10)
    days = (11, 12, 13)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        trues(length(days)),
        (gpp = falses(length(days)), reco = trues(length(days)), et = trues(length(days))),
        (gpp = values.gpp, reco = values.reco .+ 0.05f0, et = values.et .+ 0.05f0),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    parameter_names = (:k_soil10_fast, :PRIESTLEY_TAYLOR, :soildepth_evap, :soil_infil)
    parameters = ModelParameters(Float32)
    lpjml = parameters.lpjml
    theta = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end

    result = enzyme_seasonal_soil_gradient_blockwise(
        theta,
        state_factory,
        cft1,
        parameters,
        parameter_names,
        template.climate,
        days,
        template.layer_depth,
        context;
        block_days = 2,
    )

    @test isfinite(result.primal)
    @test all(isfinite, result.gradient)
end
