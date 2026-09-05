function _nitrogen_limited_seasonal_value(
    template,
    theta,
    parameter_names,
    days,
    context,
)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_seasonal_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context;
        nitrogen_limit_vcmax = true,
    )
end

@testset "Enzyme nitrogen-limited parameter optimization" begin
    template = _nitrogen_limited_transition_fixture()
    days = (11, 12, 13)
    context = ADSeasonContext(
        trues(3),
        (gpp = trues(3), reco = trues(3), et = trues(3)),
        (
            gpp = zeros(Float32, 3),
            reco = zeros(Float32, 3),
            et = zeros(Float32, 3),
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    parameter_names = (:knstore, :gmin)
    theta = Float32[cft1.knstore, cft1.gmin]
    state_factory() = begin
        state = deepcopy(template.state)
        enzyme_prepare_daily_state!(state)
        return state, enzyme_zero_tangent(state)
    end

    result = enzyme_seasonal_gradient_blockwise(
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context;
        block_days = 1,
        nitrogen_limit_vcmax = true,
    )
    @test isfinite(result.primal)
    @test all(isfinite, result.gradient)
    @test any(!iszero, result.gradient)

    finite_difference = similar(theta)
    for index in eachindex(theta)
        step = 1.0f-3 * max(abs(theta[index]), 1.0f0)
        delta = zeros(Float32, length(theta))
        delta[index] = step
        finite_difference[index] = (
            _nitrogen_limited_seasonal_value(
                template, theta .+ delta, parameter_names, days, context,
            ) -
            _nitrogen_limited_seasonal_value(
                template, theta .- delta, parameter_names, days, context,
            )
        ) / (2.0f0 * step)
    end
    @test result.gradient ≈ finite_difference rtol = 5.0f-2 atol = 5.0f-4

    direction = result.gradient / max(norm(result.gradient), eps(Float32))
    updated_theta = theta .- 1.0f-3 .* direction
    updated_loss = _nitrogen_limited_seasonal_value(
        template, updated_theta, parameter_names, days, context,
    )
    @test updated_loss < result.forward_primal

    soil_parameter_names = (
        :k_soil10_fast,
        :PRIESTLEY_TAYLOR,
        :soildepth_evap,
        :soil_infil,
    )
    lpjml = template.global_parameters.lpjml
    soil_theta = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    soil_result = enzyme_seasonal_soil_gradient_blockwise(
        soil_theta,
        state_factory,
        cft1,
        template.global_parameters,
        soil_parameter_names,
        template.climate,
        days,
        template.layer_depth,
        context;
        block_days = 1,
        nitrogen_limit_vcmax = true,
    )
    @test isfinite(soil_result.primal)
    @test all(isfinite, soil_result.gradient)
end
