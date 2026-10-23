# Preserve the production solver's primal while allowing a smooth surrogate
# to carry derivatives through its locally inactive correction. This keeps AD
# trajectories on the ordinary production path without differentiating the
# discrete bisection decisions used by the production lambda solver.
@noinline _enzyme_primal_correction(production, surrogate) = production - surrogate

function Enzyme.EnzymeRules.inactive(
    ::typeof(_enzyme_primal_correction),
    args...,
)
    return nothing
end

@inline _enzyme_lambda_primal(::Val{:surrogate}, production, surrogate) = surrogate
@inline _enzyme_lambda_primal(::Val{:production}, production, surrogate) =
    surrogate + _enzyme_primal_correction(production, surrogate)

@inline function _enzyme_c3_lambda_residual_slope(
    fac,
    lambda,
    vcmax,
    stress,
    b,
    co2,
    temperature,
    apar,
    daylength,
    lpjmlparams,
    photoparams,
)
    T = typeof(fac)
    stress < T(1e-2) && return -fac
    ko = lpjmlparams.ko25 * photoparams.q10ko ^ ((temperature - T(25)) * T(0.1))
    kc = lpjmlparams.kc25 * photoparams.q10kc ^ ((temperature - T(25)) * T(0.1))
    co2_fac = kc * (one(T) + photoparams.po2 / ko)
    tau = photoparams.tau25 * photoparams.q10tau ^ ((temperature - T(25)) * T(0.1))
    gammastar = photoparams.po2 / (T(2) * tau)
    internal_co2 = lambda * co2
    denominator_1 = internal_co2 + T(2) * gammastar
    denominator_2 = internal_co2 + co2_fac
    c1_slope = stress * lpjmlparams.alphac3 *
        co2 * T(3) * gammastar / (denominator_1 * denominator_1)
    c2_slope = co2 * (co2_fac + gammastar) / (denominator_2 * denominator_2)
    je_slope = c1_slope * apar * photoparams.cmass * photoparams.cq / daylength
    jc_slope = c2_slope * vcmax / T(24)
    je = stress * lpjmlparams.alphac3 *
        ((internal_co2 - gammastar) / denominator_1) *
        apar * photoparams.cmass * photoparams.cq / daylength
    jc = ((internal_co2 - gammastar) / denominator_2) * vcmax / T(24)
    total = je + jc
    radicand = max(zero(T), total * total - T(4) * lpjmlparams.theta * je * jc)
    root = sqrt(radicand)
    root_slope = root > zero(T) ?
        (total * (je_slope + jc_slope) -
         T(2) * lpjmlparams.theta * (je_slope * jc + je * jc_slope)) / root :
        zero(T)
    agd_slope = (je_slope + jc_slope - root_slope) /
        (T(2) * lpjmlparams.theta) * daylength
    adt_slope = agd_slope
    scale = (temperature + T(273.15)) / photoparams.p * T(8.314) / photoparams.cmass * T(1000)
    return -fac - (adt_slope > zero(T) ? adt_slope * scale : zero(T))
end

@inline function _enzyme_smooth_lambda_c3(
    fac,
    vcmax,
    stress,
    b,
    co2,
    temperature,
    apar,
    daylength,
    lpjmlparams,
    photoparams,
    upper_bound = typeof(fac)(0.85),
    max_iterations::Int = 60,
    constrain_to_upper_bound::Bool = false,
    lambda_primal_mode = Val(:surrogate),
)
    T = typeof(fac)
    production_lambda = Agrocosm.compute_lambda_c3_solution(
        fac,
        vcmax,
        stress,
        b,
        co2,
        temperature,
        apar,
        daylength,
        lpjmlparams,
        photoparams,
        upper_bound,
        max_iterations,
    )
    lambda = production_lambda
    for _ in 1:8
        residual = fac * (one(T) - lambda) - Agrocosm.c3_adtmm_scalar_impl(
            lambda,
            vcmax,
            stress,
            b,
            co2,
            temperature,
            apar,
            daylength,
            lpjmlparams,
            photoparams,
        )
        slope = _enzyme_c3_lambda_residual_slope(
            fac,
            lambda,
            vcmax,
            stress,
            b,
            co2,
            temperature,
            apar,
            daylength,
            lpjmlparams,
            photoparams,
        )
        updated_lambda = lambda - residual / slope
        lambda = constrain_to_upper_bound ?
            clamp(updated_lambda, zero(T), upper_bound) : updated_lambda
    end
    return _enzyme_lambda_primal(lambda_primal_mode, production_lambda, lambda)
end

# The production lambda wrapper packs active CFT scalars into a NamedTuple
# before launching its kernel. Keep the AD path scalar and explicit so Enzyme
# can carry the CFT tangent through the continuous fixed-iteration solve.
function _enzyme_solve_lambda_c3!(
    state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters,
    pet,
    temperature,
    co2,
    lpjmlparams,
    photoparams,
    constrain_to_upper_bound::Bool = false,
    pathway = Val(:C3),
    lambda_primal_mode = Val(:surrogate),
)
    T = eltype(Agrocosm.crop_photosynthesis_auxiliary(state).lambda)
    lambda = Agrocosm.crop_photosynthesis_auxiliary(state).lambda
    vcmax = Agrocosm.crop_photosynthesis_auxiliary(state).vcmax
    temperature_stress = Agrocosm.crop_photosynthesis_auxiliary(state).temperature_stress
    conductance = Agrocosm.crop_canopy_auxiliary(state).canopy_conductance
    fpar = Agrocosm.crop_canopy_auxiliary(state).fpar
    apar = Agrocosm.crop_canopy_auxiliary(state).apar
    for cell in eachindex(lambda)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        gpd, fac = Agrocosm.compute_canopy_water_supply(
            pet.daylength[cell], conductance[cell], T(cft.gmin), fpar[cell], co2_cell,
        )
        lambda[cell] = if gpd > T(1e-5) &&
                          temperature_stress[cell] >= T(1e-2) &&
                          pet.daylength[cell] > zero(T) && co2_cell > zero(T)
            _enzyme_weather_lambda(pathway,
                fac,
                vcmax[cell],
                temperature_stress[cell],
                T(cft.b),
                co2_cell,
                temperature[cell],
                apar[cell],
                pet.daylength[cell],
                lpjmlparams,
                photoparams,
                T(0.85),
                60,
                constrain_to_upper_bound,
                lambda_primal_mode,
            )
        else
            zero(T)
        end
    end
    return nothing
end

function _enzyme_recouple_nitrogen_water_c3!(
    state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters,
    pet,
    temperature,
    co2,
    lpjmlparams,
    photoparams,
    pathway = Val(:C3),
    lambda_primal_mode = Val(:surrogate),
)
    T = eltype(Agrocosm.crop_photosynthesis_auxiliary(state).lambda)
    photosynthesis = Agrocosm.crop_photosynthesis_auxiliary(state)
    canopy = Agrocosm.crop_canopy_auxiliary(state)
    crop = Agrocosm.crop_prognostic(state)
    assimilation = Agrocosm.crop_fluxes(state).carbon.water_limited_assimilation
    root_distribution = Agrocosm.crop_root_input(state).distribution
    soil_water = Agrocosm.soil_water_auxiliary(state).relative_content
    for cell in eachindex(photosynthesis.lambda)
        if crop.phenology.is_growing[cell] == one(eltype(crop.phenology.is_growing)) &&
           photosynthesis.lambda[cell] > zero(T) &&
           photosynthesis.temperature_stress[cell] >= T(1e-2)
            co2_cell = co2[length(co2) == 1 ? 1 : cell]
            previous_lambda = photosynthesis.lambda[cell]
            potential_conductance = canopy.canopy_conductance[cell]
            limited_conductance = Agrocosm.compute_canopy_conductance(
                assimilation[cell], co2_cell, pet.daylength[cell], canopy.fpar[cell],
                T(cft.gmin), previous_lambda,
            )
            demand = Agrocosm.compute_transpiration_demand(
                canopy.canopy_wet[cell], pet.eeq[cell], T(lpjmlparams.ALPHAM),
                T(lpjmlparams.GM), limited_conductance,
            )
            root_water = zero(T)
            for layer in axes(soil_water, 1)
                root_water += soil_water[layer, cell] * root_distribution[layer]
            end
            supply = Agrocosm.compute_transpiration_supply(
                T(cft.emax), root_water, crop.carbon.root[cell],
            ) * T(cft.fpc)
            canopy.canopy_conductance[cell] = limited_conductance

            if Agrocosm.nitrogen_water_recoupling_required(
                potential_conductance, limited_conductance, demand, supply,
            )
                constrained_conductance = Agrocosm.compute_actual_canopy_conductance(
                    limited_conductance, supply, demand, canopy.canopy_wet[cell],
                    pet.eeq[cell], T(lpjmlparams.ALPHAM), T(lpjmlparams.GM),
                )
                gpd, fac = Agrocosm.compute_canopy_water_supply(
                    pet.daylength[cell], constrained_conductance, T(cft.gmin),
                    canopy.fpar[cell], co2_cell,
                )
                if gpd > T(1e-5) && pet.daylength[cell] > zero(T) && co2_cell > zero(T)
                    photosynthesis.lambda[cell] = _enzyme_weather_lambda(pathway,
                        fac,
                        photosynthesis.vcmax[cell],
                        photosynthesis.temperature_stress[cell],
                        T(cft.b),
                        co2_cell,
                        temperature[cell],
                        canopy.apar[cell],
                        pet.daylength[cell],
                        lpjmlparams,
                        photoparams,
                        previous_lambda,
                        20,
                        true,
                        lambda_primal_mode,
                    )
                end
            end
        end
    end
    return nothing
end

function _enzyme_refresh_potential_canopy_conductance!(state, cft, pet, co2)
    T = eltype(Agrocosm.crop_photosynthesis_auxiliary(state).lambda)
    photosynthesis = Agrocosm.crop_photosynthesis_auxiliary(state)
    canopy = Agrocosm.crop_canopy_auxiliary(state)
    crop = Agrocosm.crop_prognostic(state)
    assimilation = Agrocosm.crop_fluxes(state).carbon.water_limited_assimilation
    for cell in eachindex(photosynthesis.lambda)
        if crop.phenology.is_growing[cell] == one(eltype(crop.phenology.is_growing)) &&
           photosynthesis.lambda[cell] > zero(T) &&
           photosynthesis.temperature_stress[cell] >= T(1e-2)
            co2_cell = co2[length(co2) == 1 ? 1 : cell]
            canopy.canopy_conductance[cell] = Agrocosm.compute_canopy_conductance(
                assimilation[cell], co2_cell, pet.daylength[cell], canopy.fpar[cell],
                T(cft.gmin), photosynthesis.lambda[cell],
            )
        end
    end
    return nothing
end

function _enzyme_prepare_prephenology_canopy_conductance!(
    state, cft, pet, co2, lpjmlparams,
)
    T = eltype(Agrocosm.crop_canopy_auxiliary(state).canopy_conductance)
    canopy = Agrocosm.crop_canopy_auxiliary(state)
    assimilation = Agrocosm.crop_fluxes(state).carbon.water_limited_assimilation
    for cell in eachindex(canopy.canopy_conductance)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        canopy.canopy_conductance[cell] = Agrocosm.compute_canopy_conductance(
            assimilation[cell], co2_cell, pet.daylength[cell], canopy.fpar[cell],
            T(cft.gmin), T(lpjmlparams.LAMBDA_OPT),
        )
    end
    return nothing
end

function _enzyme_finalize_nitrogen_limited_transpiration!(
    state::Agrocosm.ModelState,
    cft::Agrocosm.CFTParameters,
    pet,
    co2,
    lpjmlparams,
)
    T = eltype(Agrocosm.crop_photosynthesis_auxiliary(state).lambda)
    photosynthesis = Agrocosm.crop_photosynthesis_auxiliary(state)
    canopy = Agrocosm.crop_canopy_auxiliary(state)
    crop = Agrocosm.crop_prognostic(state)
    assimilation = Agrocosm.crop_fluxes(state).carbon.water_limited_assimilation
    transpiration = Agrocosm.crop_fluxes(state).water.transpiration_layer
    root_distribution = Agrocosm.crop_root_input(state).distribution
    soil_water = Agrocosm.soil_water_auxiliary(state)
    for cell in eachindex(photosynthesis.lambda)
        if crop.phenology.is_growing[cell] == one(eltype(crop.phenology.is_growing)) &&
           photosynthesis.lambda[cell] > zero(T) &&
           photosynthesis.temperature_stress[cell] >= T(1e-2)
            co2_cell = co2[length(co2) == 1 ? 1 : cell]
            final_conductance = Agrocosm.compute_canopy_conductance(
                assimilation[cell], co2_cell, pet.daylength[cell], canopy.fpar[cell],
                T(cft.gmin), photosynthesis.lambda[cell],
            )
            demand = Agrocosm.compute_transpiration_demand(
                canopy.canopy_wet[cell], pet.eeq[cell], T(lpjmlparams.ALPHAM),
                T(lpjmlparams.GM), final_conductance,
            )
            root_water = zero(T)
            for layer in axes(soil_water.relative_content, 1)
                root_water += soil_water.relative_content[layer, cell] *
                    root_distribution[layer]
            end
            transpiration_scale = root_water > zero(T) ?
                demand * canopy.fpar[cell] / root_water : zero(T)
            for layer in axes(soil_water.relative_content, 1)
                transpiration[layer, cell], _ = Agrocosm.compute_layer_transpiration(
                    transpiration_scale, root_distribution[layer],
                    soil_water.relative_content[layer, cell],
                    soil_water.holding_capacity_storage[layer, cell],
                )
            end
            canopy.canopy_conductance[cell] = final_conductance
        end
    end
    return nothing
end
