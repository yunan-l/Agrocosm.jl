using LinearAlgebra

const _AgrocosmEnzymeExt = Base.get_extension(Agrocosm, :AgrocosmEnzymeExt)

function _management_split_event_mineral_n(
    theta, state, context, day, lpjmlparams,
)
    _AgrocosmEnzymeExt._apply_management_split_fertilizer!(
        state, theta[1], theta[2], context, day, lpjmlparams,
    )
    nitrogen = Agrocosm.soil_nitrogen_prognostic(state)
    return nitrogen.nitrate[1, 1] + nitrogen.ammonium[1, 1]
end

function _management_adaptation_state_factory(climate_days::Int = 18)
    template = _daily_transition_fixture(climate_days, 0)
    state = deepcopy(template.state)
    daily_crop_C3!(
        1, 10, ProcessModules(cft1, template.global_parameters),
        template.climate, state;
        fertilizer = :no,
        manure = false,
        with_tillage = true,
        update_vernalization_requirement = false,
        reuse_output = true,
    )
    nitrogen = Agrocosm.crop_prognostic(state).nitrogen
    fill!(nitrogen.pending_fertilizer, 0.0f0)
    fill!(nitrogen.pending_manure, 0.0f0)
    enzyme_prepare_daily_state!(state)
    return merge(template, (; state, day = 11))
end

function _management_adaptation_loss(template, theta, context; kwargs...)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_management_yield_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context;
        kwargs...,
    )
end

function _management_split_adaptation_loss(template, theta, context, days = 11:18)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_management_yield_split_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        days,
        template.layer_depth,
        context,
    )
end

function _management_joint_adaptation_loss(template, theta, context, days = 11:18)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_joint_adaptation_yield_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        days,
        template.layer_depth,
        context,
    )
end

@testset "Enzyme fixed-event management adaptation" begin
    template = _management_adaptation_state_factory()
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55]
    direction = Float32[1]
    crop_nitrogen = Agrocosm.crop_prognostic(template.state).nitrogen
    crop_nitrogen_flux = Agrocosm.crop_fluxes(template.state).nitrogen
    @test iszero(only(crop_nitrogen.pending_fertilizer))
    @test iszero(only(crop_nitrogen.pending_manure))
    @test iszero(only(crop_nitrogen_flux.prescribed_fertilizer_input))
    @test iszero(only(crop_nitrogen_flux.prescribed_manure_input))
    primal = _management_adaptation_loss(template, theta, context)
    @test isfinite(primal)
    @test primal == _management_adaptation_loss(
        template, theta, context; crop_resp_fix = true,
    )
    dynamic_respiration_primal = _management_adaptation_loss(
        template, theta, context; crop_resp_fix = false,
    )
    @test isfinite(dynamic_respiration_primal)

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    ad = enzyme_forward_directional(
        enzyme_management_yield_loss,
        theta,
        direction,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context;
        return_primal = true,
    )
    @test ad.primal ≈ primal rtol = 1.0f-6 atol = 1.0f-7

    step = 1.0f-3
    finite_difference = (
        _management_adaptation_loss(template, theta .+ step .* direction, context) -
        _management_adaptation_loss(template, theta .- step .* direction, context)
    ) / (2.0f0 * step)
    @test ad.directional ≈ finite_difference rtol = 5.0f-2 atol = 1.0f-4

    state_with_latent_dose = deepcopy(template.state)
    state_without_latent_dose = deepcopy(template.state)
    pending = Agrocosm.crop_prognostic(state_with_latent_dose).nitrogen.pending_fertilizer
    pending_amount = 7.5f0
    fill!(pending, pending_amount)
    fill!(Agrocosm.crop_prognostic(state_without_latent_dose).nitrogen.pending_fertilizer, 0.0f0)
    enzyme_prepare_daily_state!(state_with_latent_dose)
    enzyme_prepare_daily_state!(state_without_latent_dose)
    loss_with_latent_dose = enzyme_management_yield_loss(
        theta, state_with_latent_dose, cft1, template.global_parameters,
        template.climate, 11:18, template.layer_depth, context;
        nitrogen_limit_vcmax = true,
    )
    loss_without_latent_dose = enzyme_management_yield_loss(
        theta, state_without_latent_dose, cft1, template.global_parameters,
        template.climate, 11:18, template.layer_depth, context;
        nitrogen_limit_vcmax = true,
    )
    @test loss_with_latent_dose ≈ loss_without_latent_dose rtol = 1.0f-6
    @test only(pending) == pending_amount

    @test_throws ArgumentError ManagementAdaptationContext((11, 11), (0.2f0, 0.8f0))
    @test_throws ArgumentError ManagementAdaptationContext((11, 14), (0.1f0, 0.8f0))
end


@testset "Enzyme joint management and cultivar adaptation" begin
    template = _management_adaptation_state_factory(26)
    days = 11:26
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55, 0.2, 1.0, 0.0, 1.0, 1.0]

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    gradient = enzyme_forward_gradient(
        enzyme_joint_adaptation_yield_loss,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        days,
        template.layer_depth,
        context,
    )
    @test all(isfinite, gradient)
    @test abs(gradient[2]) > 1.0f-3

    storage_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(storage_state)
    enzyme_joint_adaptation_yield_loss(
        theta, storage_state, cft1, template.global_parameters,
        template.climate, days, template.layer_depth, context,
    )
    @test only(Agrocosm.crop_prognostic(storage_state).carbon.storage) > 1.0f-3

    steps = Float32[1.0f-2, 1.0f-2, 1.0f-3, 1.0f-1, 1.0f-3, 1.0f-3]
    finite_difference = similar(theta)
    for index in eachindex(theta)
        perturbation = zeros(Float32, length(theta))
        perturbation[index] = steps[index]
        finite_difference[index] = (
            _management_joint_adaptation_loss(
                template, theta .+ perturbation, context, days,
            ) -
            _management_joint_adaptation_loss(
                template, theta .- perturbation, context, days,
            )
        ) / (2.0f0 * steps[index])
    end
    @test all(isfinite, finite_difference)
    @test gradient ≈ finite_difference rtol = 5.0f-2 atol = 2.0f-4
end


@testset "Enzyme two-event fertilizer pulse tangent" begin
    template = _management_adaptation_state_factory()
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0))
    theta = Float32[24.55, 0.2]
    direction = Float32[0, 1]

    for (day, expected) in ((11, theta[1]), (14, -theta[1]), (12, 0.0f0))
        state_factory() = begin
            state = deepcopy(template.state)
            return state, enzyme_zero_tangent(state)
        end
        directional = enzyme_forward_directional(
            _management_split_event_mineral_n,
            theta,
            direction,
            state_factory,
            context,
            day,
            template.global_parameters.lpjml,
        )
        @test directional ≈ expected rtol = 1.0f-6 atol = 1.0f-6
    end
end


@testset "Enzyme two-event fertilizer allocation" begin
    template = _management_adaptation_state_factory(26)
    days = 11:26
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55, 0.2]

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    gradient = enzyme_forward_gradient(
        enzyme_management_yield_split_loss,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        days,
        template.layer_depth,
        context,
    )
    @test all(isfinite, gradient)
    @test abs(gradient[2]) > 1.0f-3

    storage_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(storage_state)
    enzyme_management_yield_split_loss(
        theta, storage_state, cft1, template.global_parameters,
        template.climate, days, template.layer_depth, context,
    )
    @test only(Agrocosm.crop_prognostic(storage_state).carbon.storage) > 1.0f-3

    for index in eachindex(theta)
        step = 1.0f-2
        direction = zeros(Float32, length(theta))
        direction[index] = step
        finite_difference = (
            _management_split_adaptation_loss(
                template, theta .+ direction, context, days,
            ) -
            _management_split_adaptation_loss(
                template, theta .- direction, context, days,
            )
        ) / (2.0f0 * step)
        @test gradient[index] ≈ finite_difference rtol = 5.0f-2 atol = 2.0f-4
    end

    three_event_context = ManagementAdaptationContext(
        (11, 13, 15), (0.2f0, 0.3f0, 0.5f0),
    )
    @test_throws ArgumentError _management_split_adaptation_loss(
        template, theta, three_event_context,
    )
end
