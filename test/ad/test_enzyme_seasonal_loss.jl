function _seasonal_production_values(template, days)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    processes = ProcessModules(template.base_cft, template.global_parameters)
    gpp = zeros(Float32, length(days))
    reco = zeros(Float32, length(days))
    et = zeros(Float32, length(days))
    daily_crop! = template.base_cft.path == 1 ? daily_crop_C3! : daily_crop_C4!
    for index in eachindex(days)
        daily_crop!(
            days[index], days[index], processes,
            template.climate, state;
            fertilizer = :yes,
            manure = true,
            with_tillage = true,
            update_vernalization_requirement = false,
            reuse_output = true,
        )
        values = _production_daily_transition_value(state, :gpp)
        gpp[index] = values
        reco[index] = _production_daily_transition_value(state, :reco)
        et[index] = _production_daily_transition_value(state, :et)
    end
    return (; gpp, reco, et)
end

function _seasonal_loss_value(template, theta, parameter_names, days, context)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_seasonal_loss(
        theta,
        state,
        template.base_cft,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context,
    )
end

function _expected_seasonal_loss(values, context)
    losses = zeros(Float32, 3)
    names = (:gpp, :reco, :et)
    for index in eachindex(context.growth_mask)
        context.growth_mask[index] || continue
        for target_index in eachindex(names)
            target = names[target_index]
            getfield(context.valid_masks, target)[index] || continue
            residual = (
                getfield(values, target)[index] - getfield(context.observations, target)[index]
            ) / getfield(context.scales, target)
            losses[target_index] += residual * residual
        end
    end
    return losses[1] / max(context.counts.gpp, 1) +
        losses[2] / max(context.counts.reco, 1) +
        losses[3] / max(context.counts.et, 1)
end

@testset "Enzyme seasonal loss" begin
    template = _daily_transition_fixture()
    days = (11, 12, 13)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        Bool[false, true, true],
        (gpp = trues(3), reco = trues(3), et = trues(3)),
        (
            gpp = values.gpp .+ Float32[999, 1, -1],
            reco = values.reco .+ Float32[999, 0.5, -0.5],
            et = values.et .+ Float32[999, 0.5, -0.5],
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)
    expected = _expected_seasonal_loss(values, context)
    adapter_value = _seasonal_loss_value(template, theta, parameter_names, days, context)
    @test adapter_value ≈ expected rtol = 1.0f-2 atol = 5.0f-3

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    ad = enzyme_forward_directional(
        enzyme_seasonal_loss,
        theta,
        direction,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context;
        return_primal = true,
    )
    @test ad.primal ≈ adapter_value rtol = 1.0f-6 atol = 1.0f-7

    gradient = enzyme_forward_gradient(
        enzyme_seasonal_loss,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context,
    )
    projection = dot(gradient, direction)
    projection_scale = sum(abs, gradient .* direction)
    # Float32 coordinate and directional passes accumulate in different orders.
    @test abs(ad.directional - projection) <=
        5.0f-5 * projection_scale + 1.0f-6

    finite_difference = similar(theta)
    for index in eachindex(theta)
        step = 1.0f-3 * max(abs(theta[index]), 1.0f0)
        delta = zeros(Float32, length(theta))
        delta[index] = step
        finite_difference[index] = (
            _seasonal_loss_value(template, theta .+ delta, parameter_names, days, context) -
            _seasonal_loss_value(template, theta .- delta, parameter_names, days, context)
        ) / (2.0f0 * step)
    end
    @test gradient ≈ finite_difference rtol = 2.0f-2 atol = 1.0f-5
end

@testset "Enzyme C4 seasonal pathway" begin
    template = _daily_transition_fixture(; base_cft = cft3)
    days = (template.day,)
    values = _seasonal_production_values(template, days)
    @test values.gpp[1] > 0.0f0
    context = ADSeasonContext(
        trues(1),
        (gpp = trues(1), reco = trues(1), et = trues(1)),
        (
            gpp = values.gpp .+ 1.0f0,
            reco = values.reco .+ 0.5f0,
            et = values.et .+ 0.5f0,
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    cft_names = (:gmin, :b)
    cft_theta = Float32[cft3.gmin, cft3.b]
    soil_names = (:soildepth_evap,)
    soil_theta = Float32[template.global_parameters.lpjml.soildepth_evap]
    expected = _expected_seasonal_loss(values, context)
    cft_loss = _seasonal_loss_value(template, cft_theta, cft_names, days, context)
    @test cft_loss ≈ expected rtol = 1.0f-2 atol = 5.0f-3

    daily_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(daily_state)
    daily_gpp = enzyme_daily_transition_objective(
        cft_theta, daily_state, cft3, template.global_parameters,
        template.climate, cft_names, template.day, :gpp, template.layer_depth,
    )
    @test daily_gpp ≈ values.gpp[1] rtol = 5.0f-5 atol = 1.5f-4

    soil_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(soil_state)
    soil_loss = enzyme_seasonal_soil_loss(
        soil_theta, soil_state, cft3, template.global_parameters, soil_names,
        template.climate, days, template.layer_depth, context,
    )
    @test soil_loss ≈ cft_loss rtol = 1.0f-5 atol = 1.0f-6

    joint_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(joint_state)
    joint_loss = enzyme_seasonal_joint_loss(
        cft_theta, soil_theta, joint_state, cft3, template.global_parameters,
        cft_names, soil_names, template.climate, days, template.layer_depth, context,
    )
    @test joint_loss ≈ cft_loss rtol = 1.0f-5 atol = 1.0f-6

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    blockwise = enzyme_seasonal_gradient_blockwise(
        cft_theta, state_factory, cft3, template.global_parameters,
        template.climate, cft_names, days, template.layer_depth, context;
        block_days = 1,
    )
    @test blockwise.forward_primal ≈ cft_loss rtol = 1.0f-6 atol = 1.0f-6
    @test blockwise.primal ≈ blockwise.forward_primal rtol = 2.0f-5 atol = 2.0f-6
    @test all(isfinite, blockwise.gradient)
end
