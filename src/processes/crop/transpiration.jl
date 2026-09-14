"""
transpiration!(photos_adtmm, CFT, crop, pet, soil, co2; lpjmlparams=lpjmlparams)

Compute water demand/supply balance and layer-resolved transpiration uptake.
"""
function transpiration!(photos_adtmm::AbstractArray{T},
                        CFT::CFTParameters,
                        crop,
                        pet::PetPar,
                        soil,
                        co2::AbstractArray{T};
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        use_precomputed_conductance::Bool = false,
) where {T <: AbstractFloat}

    # Root-zone weighted water availability is accumulated inside the cell
    # kernel, avoiding a separate broadcast and reduction array every day.
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop_prognostic(crop).carbon.root))
    # demand = ifelse.(crop_canopy_auxiliary(crop).canopy_conductance .> 0, (1 .- crop_canopy_auxiliary(crop).canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop_canopy_auxiliary(crop).canopy_conductance), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    kernel_params = (
        lpjmlparams = lpjmlparams,
        soil_layers = 5,
        use_precomputed_conductance = use_precomputed_conductance,
        canopy_height = T(fao56_canopy_height(CFT.name)),
    )
    weather = weather_input(crop)

    launch_1D!(water_demand_supply_kernel!,
               crop_canopy_auxiliary(crop).canopy_conductance,
               photos_adtmm,
               co2,
               pet.daylength,
               crop_canopy_auxiliary(crop).fpar,
               crop_fluxes(crop).water.transpiration_layer,
               crop_prognostic(crop).water.demand_sum,
               crop_prognostic(crop).water.supply_sum,
               crop_stress_auxiliary(crop).water_deficit,
               crop_prognostic(crop).water.sufficiency,
               crop_prognostic(crop).carbon.root,
               crop_canopy_auxiliary(crop).canopy_wet,
               crop_prognostic(crop).phenology.is_growing,
               pet.eeq,
               weather.temp,
               weather.wind,
               weather.vapour_deficit,
               crop_root_input(crop).distribution,
               crop_root_auxiliary(crop).zone_available_water,
               soil_water_auxiliary(soil).relative_content,
               soil_water_auxiliary(soil).holding_capacity_storage,
               CFT,
               kernel_params)

end

"""
    prepare_prephenology_canopy_conductance!(CFT, crop, daylength, co2)

Store LPJmL's raw `gp_sum` crop conductance after photosynthesis has been
evaluated with the canopy state present before today's phenology update.
"""
function prepare_prephenology_canopy_conductance!(
    CFT::CFTParameters,
    crop,
    daylength::AbstractArray{T},
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    launch_1D!(
        prepare_prephenology_canopy_conductance_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        co2,
        daylength,
        crop_canopy_auxiliary(crop).fpar,
        T(CFT.gmin),
        T(lpjmlparams.LAMBDA_OPT),
    )
    return nothing
end

@kernel inbounds = true function prepare_prephenology_canopy_conductance_kernel!(
    conductance::AbstractArray{T},
    assimilation::AbstractArray{T},
    co2::AbstractArray{T},
    daylength::AbstractArray{T},
    fpar::AbstractArray{T},
    gmin::T,
    lambda_optimum::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    co2_cell = co2[length(co2) == 1 ? 1 : cell]
    conductance[cell] = compute_canopy_conductance(
        assimilation[cell], co2_cell, daylength[cell], fpar[cell], gmin,
        lambda_optimum,
    )
end

"""Return whether LPJmL repeats the water-limited lambda solve after N stress."""
@inline nitrogen_water_recoupling_required(
    potential_conductance::T,
    nitrogen_limited_conductance::T,
    demand::T,
    supply::T,
) where {T <: AbstractFloat} =
    potential_conductance - nitrogen_limited_conductance > T(0.01) &&
    demand - supply > T(0.1)

"""Store LPJmL's pre-N-limitation `gc_new` from potential assimilation."""
function refresh_potential_canopy_conductance!(
    CFT::CFTParameters,
    crop,
    daylength::AbstractArray{T},
    co2::AbstractArray{T},
) where {T <: AbstractFloat}
    launch_1D!(
        refresh_potential_canopy_conductance_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).fpar,
        crop_prognostic(crop).phenology.is_growing,
        daylength,
        co2,
        T(CFT.gmin),
    )
    return nothing
end

@kernel inbounds = true function refresh_potential_canopy_conductance_kernel!(
    conductance::AbstractArray{T},
    potential_assimilation::AbstractArray{T},
    lambda::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    fpar::AbstractArray{T},
    is_growing::AbstractArray{S},
    daylength::AbstractArray{T},
    co2::AbstractArray{T},
    gmin::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        conductance[cell] = compute_canopy_conductance(
            potential_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, lambda[cell],
        )
    end
end

"""
    recouple_nitrogen_water!(pathway, CFT, crop, pet, soil, temperature, co2)

Apply LPJmL's conditional second lambda solve after leaf nitrogen has limited
Vcmax. The first water pass remains responsible for seasonal water-stress
diagnostics; this correction only updates today's conductance and lambda.
"""
function recouple_nitrogen_water!(
    pathway::Union{Val{:C3}, Val{:C4}},
    CFT::CFTParameters,
    crop,
    pet::PetPar,
    soil,
    temperature::AbstractArray{T},
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
    photoparams::PhotoParams = photoparams,
) where {T <: AbstractFloat}
    kernel_params = (
        pathway = pathway,
        b = T(CFT.b),
        emax = T(CFT.emax),
        depletion_fraction = T(CFT.depletion_fraction),
        depletion_demand_slope = T(CFT.depletion_demand_slope),
        fpc = T(CFT.fpc),
        gmin = T(CFT.gmin),
        soil_layers = 5,
        canopy_height = T(fao56_canopy_height(CFT.name)),
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )
    weather = weather_input(crop)
    launch_1D!(
        nitrogen_water_recoupling_kernel!,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).apar,
        crop_canopy_auxiliary(crop).canopy_wet,
        crop_prognostic(crop).phenology.is_growing,
        crop_prognostic(crop).carbon.root,
        crop_root_input(crop).distribution,
        soil_water_auxiliary(soil).relative_content,
        pet.daylength,
        pet.eeq,
        temperature,
        weather.wind,
        weather.vapour_deficit,
        co2,
        kernel_params,
    )
    return nothing
end

@kernel inbounds = true function nitrogen_water_recoupling_kernel!(
    lambda::AbstractArray{T},
    vcmax::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    limited_assimilation::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    canopy_wet::AbstractArray{T},
    is_growing::AbstractArray{S},
    root_carbon::AbstractArray{T},
    root_distribution::AbstractArray{T},
    soil_water::AbstractArray{M},
    daylength::AbstractArray{T},
    equilibrium_evaporation::AbstractArray{T},
    temperature::AbstractArray{T},
    wind::AbstractArray{T},
    vapour_deficit::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    @unpack pathway, b, emax, fpc, gmin, soil_layers, lpjmlparams, photoparams = kernel_params
    @unpack depletion_fraction, depletion_demand_slope = kernel_params
    @unpack ALPHAM, GM = lpjmlparams
    canopy_height = T(kernel_params.canopy_height)
    # Penman-Monteith's stomatal sensitivity, expressed as LPJmL's own two
    # constants so that no equation downstream changes. `aerodynamic_coupling`
    # is 0 by default, which returns those constants untouched.
    alpha_c, shape_c = coupled_demand_parameters(
        T(ALPHAM), T(GM), equilibrium_evaporation[cell], temperature[cell],
        vapour_deficit[cell], wind[cell], canopy_height * fpar[cell],
        T(lpjmlparams.aerodynamic_coupling),
    )

    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        previous_lambda = lambda[cell]
        potential_conductance = conductance[cell]
        limited_conductance = compute_canopy_conductance(
            limited_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, previous_lambda,
        )
        demand = compute_transpiration_demand(
            canopy_wet[cell], equilibrium_evaporation[cell], alpha_c, shape_c,
            limited_conductance,
        )
        root_water = zero(T)
        for layer in 1:soil_layers
            root_water += soil_water[layer, cell] * root_distribution[layer]
        end
        supply = compute_transpiration_supply(
            emax, root_water, root_carbon[cell],
            demand_adjusted_depletion(
                depletion_fraction,
                equilibrium_evaporation[cell] * alpha_c,
                depletion_demand_slope,
            ),
        ) * fpc
        conductance[cell] = limited_conductance

        if nitrogen_water_recoupling_required(
            potential_conductance, limited_conductance, demand, supply,
        )
            constrained_conductance = compute_actual_canopy_conductance(
                limited_conductance, supply, demand, canopy_wet[cell],
                equilibrium_evaporation[cell], alpha_c, shape_c,
            )
            gpd, fac = compute_canopy_water_supply(
                daylength[cell], constrained_conductance, gmin, fpar[cell], co2_cell,
            )
            if gpd > T(1e-5) && daylength[cell] > zero(T) && co2_cell > zero(T)
                if pathway isa Val{:C3}
                    lambda[cell] = compute_lambda_c3_solution(
                        fac, vcmax[cell], temperature_stress[cell], b, co2_cell,
                        temperature[cell], apar[cell], daylength[cell], lpjmlparams,
                        photoparams, previous_lambda, 20,
                    )
                else
                    lambda[cell] = compute_lambda_c4_solution(
                        fac, vcmax[cell], temperature_stress[cell], b,
                        temperature[cell], apar[cell], daylength[cell], lpjmlparams,
                        photoparams, previous_lambda, 20,
                    )
                end
            end
        end
    end
end

"""Overwrite today's transpiration layers using final N-limited photosynthesis."""
function finalize_nitrogen_limited_transpiration!(
    CFT::CFTParameters,
    crop,
    pet::PetPar,
    soil,
    co2::AbstractArray{T};
    lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    kernel_params = (gmin = T(CFT.gmin), soil_layers = 5,
                     canopy_height = T(fao56_canopy_height(CFT.name)),
                     lpjmlparams = lpjmlparams)
    weather = weather_input(crop)
    launch_1D!(
        finalize_nitrogen_limited_transpiration_kernel!,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_fluxes(crop).water.transpiration_layer,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).canopy_wet,
        crop_prognostic(crop).phenology.is_growing,
        crop_root_input(crop).distribution,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).holding_capacity_storage,
        pet.daylength,
        pet.eeq,
        weather.temp,
        weather.wind,
        weather.vapour_deficit,
        co2,
        kernel_params,
    )
    return nothing
end

@kernel inbounds = true function finalize_nitrogen_limited_transpiration_kernel!(
    conductance::AbstractArray{T},
    transpiration_layer::AbstractArray{T},
    limited_assimilation::AbstractArray{T},
    lambda::AbstractArray{T},
    temperature_stress::AbstractArray{T},
    fpar::AbstractArray{T},
    canopy_wet::AbstractArray{T},
    is_growing::AbstractArray{S},
    root_distribution::AbstractArray{T},
    soil_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    daylength::AbstractArray{T},
    equilibrium_evaporation::AbstractArray{T},
    air_temperature::AbstractArray{T},
    wind::AbstractArray{T},
    vapour_deficit::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    @unpack gmin, soil_layers, lpjmlparams = kernel_params
    @unpack ALPHAM, GM = lpjmlparams
    canopy_height = T(kernel_params.canopy_height)
    # Penman-Monteith's stomatal sensitivity, expressed as LPJmL's own two
    # constants so that no equation downstream changes. `aerodynamic_coupling`
    # is 0 by default, which returns those constants untouched.
    alpha_c, shape_c = coupled_demand_parameters(
        T(ALPHAM), T(GM), equilibrium_evaporation[cell], air_temperature[cell],
        vapour_deficit[cell], wind[cell], canopy_height * fpar[cell],
        T(lpjmlparams.aerodynamic_coupling),
    )

    if is_growing[cell] == one(S) && lambda[cell] > zero(T) &&
       temperature_stress[cell] >= T(1e-2)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        final_conductance = compute_canopy_conductance(
            limited_assimilation[cell], co2_cell, daylength[cell], fpar[cell],
            gmin, lambda[cell],
        )
        demand = compute_transpiration_demand(
            canopy_wet[cell], equilibrium_evaporation[cell], alpha_c, shape_c,
            final_conductance,
        )
        root_water = zero(T)
        for layer in 1:soil_layers
            root_water += soil_water[layer, cell] * root_distribution[layer]
        end
        transpiration = root_water > zero(T) ? demand * fpar[cell] / root_water : zero(T)
        for layer in 1:soil_layers
            transpiration_layer[layer, cell], _ = compute_layer_transpiration(
                transpiration, root_distribution[layer], soil_water[layer, cell],
                holding_storage[layer, cell],
            )
        end
        conductance[cell] = final_conductance
    end
end

"""Compute LPJmL canopy conductance from water-limited assimilation."""
@inline function compute_canopy_conductance(
    water_limited_assimilation::T,
    co2::T,
    daylength::T,
    fpar::T,
    minimum_conductance::T,
    lambda_optimum::T,
) where {T <: AbstractFloat}
    co2_bar = co2 * T(1e-5)
    co2_bar > zero(T) && daylength > zero(T) || return zero(T)
    denominator = co2_bar * (one(T) - lambda_optimum) * hour2sec(daylength)
    return T(1.6) * water_limited_assimilation / denominator +
           minimum_conductance * fpar
end

"""
    compute_available_fraction(root_water, depletion_fraction)

Root-weighted relative soil water with the readily-available plateau applied.

LPJmL's transpiration supply is LINEAR in `root_water`, the root-weighted relative
plant-available water content, which equals 1 only at field capacity. Stress
therefore begins the moment the soil starts drying and never stops: measured
globally, mean season water sufficiency runs 75.7 to 91.1 across the entire
precipitation distribution and never reaches 100 even in the wettest bin, so
every additional millimetre still buys yield. That is the source of the model's
excess yield variance - 2.9 to 3.5 times the US county statistics - and of its
over-determination by rainfall. See `docs/22`.

The linear form is the land-surface convention. The crop-model convention is
FAO-56's readily available water, `RAW = p * TAW`, with no stress until depletion
exceeds it, so the first `p` of available water is extracted at the potential
rate. APSIM, DSSAT, CropSyst, EPIC and SWAP all carry a form of that threshold.

`depletion_fraction = 0` reproduces the LPJmL form BITWISE - `min(1, wr / 1)` is
`wr` for `wr <= 1`, and `relative_water` is clamped to [0, 1] at source
(`water_ice_pools.jl:36`) - which is the ablation contract every mechanism in
this project ships with.
"""
@inline function compute_available_fraction(
    root_water::T, depletion_fraction::T,
) where {T <: AbstractFloat}
    plateau = one(T) - depletion_fraction
    # A degenerate `p = 1` would mean "never stressed"; division would give Inf.
    plateau > zero(T) || return root_water > zero(T) ? one(T) : zero(T)
    return min(one(T), root_water / plateau)
end

"""
    demand_adjusted_depletion(depletion_fraction, potential_et, slope)

FAO-56's own adjustment of `p` for evaporative demand.

Table 22's values apply at a crop evapotranspiration of about 5 mm/day, and the
note beneath it gives `p_adjusted = p + 0.04 * (5 - ET_c)`, bounded to
[0.1, 0.8]. A crop under a thirsty atmosphere reaches the end of its readily
available water sooner, so the plateau is SHORTER there and longer in a humid
climate. That is the same document the tabulated values come from, so it is
literature rather than a fit.

It is also the reason a single global `p` cannot generalise. `docs/29` measured
the same failure for temperature from the other side: one absolute threshold
applied everywhere puts its response where the observed sensitivity is not. Here
the correction is published rather than inferred.

`potential_et` is the Priestley-Taylor POTENTIAL, `eeq * ALPHAM`, and not the
model's own `demand`. FAO-56's ET_c is the evapotranspiration of a WELL-WATERED
crop; the model's demand is computed through canopy conductance, so it collapses
under exactly the stress an arid cell should be reporting. Measured on five real
cells: `demand` averages 1.17 to 2.98 mm/day and the potential 2.84 to 6.04, so
`demand` sits below FAO-56's 5 mm/day reference nearly always and the adjustment
could only ever LENGTHEN the plateau - a near-uniform shift with no spatial
contrast, which is what the first global arm measured (every regional column
within 0.003 of the unadjusted arm). The potential is centred on the reference
and spans it in both directions, which is the contrast the adjustment exists to
express.

Skipped entirely when the slope is zero OR the plateau is off, rather than
evaluated with a zero slope: the [0.1, 0.8] bound would otherwise turn
`depletion_fraction = 0` into 0.1 and silently break the ablation contract that
every mechanism in this project ships with.
"""
@inline function demand_adjusted_depletion(
    depletion_fraction::T, potential_et::T, slope::T,
) where {T <: AbstractFloat}
    (slope > zero(T) && depletion_fraction > zero(T)) || return depletion_fraction
    return clamp(depletion_fraction + slope * (T(5) - potential_et), T(0.1), T(0.8))
end

"""Compute root-water-limited transpiration supply for one crop column."""
@inline compute_transpiration_supply(
    maximum_supply::T, root_water::T, root_carbon::T,
    depletion_fraction::T = zero(T),
) where {T <: AbstractFloat} =
    maximum_supply * compute_available_fraction(root_water, depletion_fraction) *
    (one(T) - exp(T(-0.04) * root_carbon))

"""
    aerodynamic_conductance(wind, canopy_height, reference_height)

Canopy aerodynamic conductance in mm s⁻¹ from the logarithmic wind profile, with
`d = 0.67h` and `z0 = 0.123h` (Monteith and Unsworth).

The height enters only inside a logarithm, so the development scaling below -
`fpar` in place of a tracked stem height, which this model does not carry - costs
about 20% in `ga` for a factor-two error in height, against the factor of 3.5 in
stomatal sensitivity this whole term exists to correct.
"""
@inline function aerodynamic_conductance(wind::T,
                                         canopy_height::T,
                                         reference_height::T) where {T <: AbstractFloat}
    height = max(canopy_height, T(0.05))
    displacement = T(0.67) * height
    roughness = T(0.123) * height
    level = max(reference_height, height + one(T))
    speed = max(wind, T(0.5))
    friction = T(0.41) * speed / log((level - displacement) / roughness)
    return T(1000) * friction * friction / speed
end

"""
    coupled_demand_parameters(alpha, conductance_shape, equilibrium_evaporation,
                              temperature, vapour_deficit, wind, canopy_height,
                              rate)

Return `(alpha, conductance_shape)` blended toward their Penman-Monteith values.

Penman-Monteith and LPJmL's demand are the SAME function of stomatal conductance.
Dividing PM through by `Δ/γ + 1` gives
`E = (N/M) * gc / (gc + ga/M)` with `N = eeq*M + K*D*ga`, `M = Δ/γ + 1`,
which is `alpha*eeq*gc/(gc + GM*alpha)` for `alpha = N/(M*eeq)` and
`GM = ga/(M*alpha)`. So the coupled demand needs no new equation anywhere - only
these two numbers, computed per cell per day, in place of two constants. That is
why the blend is exact at both ends and why `rate = 0` returns the shipped
constants unchanged.

`K = 86400*rho_cp/(lambda*gamma)`, the imposed-evaporation coefficient, is 0.6477
mm day⁻¹ per Pa per m s⁻¹ at sea level. `gamma` is held at its sea-level value:
the model carries no surface pressure, and `Δ/γ` moves by under 3% over the
elevation range of the world's cropland.
"""
@inline function coupled_demand_parameters(alpha::T,
                                           conductance_shape::T,
                                           equilibrium_evaporation::T,
                                           temperature::T,
                                           vapour_deficit::T,
                                           wind::T,
                                           canopy_height::T,
                                           rate::T) where {T <: AbstractFloat}
    blend = clamp(rate, zero(T), one(T))
    (blend > zero(T) && vapour_deficit > zero(T) && canopy_height > zero(T) &&
     equilibrium_evaporation > zero(T)) || return (alpha, conductance_shape)
    deficit = vapour_deficit
    slope = saturation_vapour_pressure_slope(temperature)
    # Pascals throughout, the convention `saturation_vapour_pressure` and the
    # P-model already use: gamma = 66.5 Pa/K and the imposed-evaporation
    # coefficient 86400*rho_cp/(lambda*gamma) = 0.6477 mm/day per Pa per m/s.
    ratio = slope / T(66.5) + one(T)
    ga = aerodynamic_conductance(wind, canopy_height, T(2)) / T(1000)
    numerator = equilibrium_evaporation * ratio + T(0.6477) * deficit * ga
    coupled_alpha = numerator / (ratio * equilibrium_evaporation)
    coupled_shape = T(1000) * ga / (ratio * coupled_alpha)
    return (alpha + blend * (coupled_alpha - alpha),
            conductance_shape + blend * (coupled_shape - conductance_shape))
end

"""Compute transpiration demand with LPJmL's 0.99 water-stress wetness cap."""
@inline function compute_transpiration_demand(
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
    conductance::T,
) where {T <: AbstractFloat}
    conductance > zero(T) || return zero(T)
    # Interception/soil evaporation retain the original stand wetness (wet_all
    # in LPJmL). Only water_stressed uses a canopy wetness capped at 0.99.
    return (one(T) - min(canopy_wet, T(0.99))) * equilibrium_evaporation * alpha /
           (one(T) + (conductance_shape * alpha) / conductance)
end

"""Return the LPJmL 0–100 seasonal water-sufficiency diagnostic."""
@inline function compute_water_sufficiency(
    supplied::T, demanded::T,
) where {T <: AbstractFloat}
    demanded > zero(T) || return T(100)
    return clamp(T(100) * supplied / demanded, zero(T), T(100))
end

"""Cap one layer's transpiration extraction by its plant-available water."""
@inline function compute_layer_transpiration(
    transpiration::T,
    root_fraction::T,
    relative_water::T,
    holding_storage::T,
) where {T <: AbstractFloat}
    unconstrained = transpiration * root_fraction * relative_water
    capacity = relative_water * holding_storage
    return min(unconstrained, capacity), unconstrained > capacity
end

"""Recover actual canopy conductance after layer-wise uptake capping."""
@inline function compute_actual_canopy_conductance(
    current_conductance::T,
    actual_supply::T,
    demand::T,
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
) where {T <: AbstractFloat}
    actual_supply < demand && equilibrium_evaporation > zero(T) || return current_conductance
    denominator = (one(T) - min(canopy_wet, T(0.99))) *
                  equilibrium_evaporation * alpha - actual_supply
    return denominator > zero(T) ? conductance_shape * alpha * actual_supply / denominator : zero(T)
end

@kernel inbounds = true function water_demand_supply_kernel!(
                                             crop_gp::AbstractArray{T},
                                             photos_adtmm::AbstractArray{T},
                                             co2::AbstractArray{T},
                                             daylength::AbstractArray{T},
                                             crop_fpar::AbstractArray{T},
                                             crop_trans_layer::AbstractArray{T},
                                             crop_w_demandsum::AbstractArray{T},
                                             crop_w_supplysum::AbstractArray{T},
                                             crop_wdf::AbstractArray{T},
                                             crop_wscal::AbstractArray{T},
                                             crop_rootc::AbstractArray{T},
                                             crop_canopy_wet::AbstractArray{T},
                                             crop_isgrowing::AbstractArray{S},
                                             pet_eeq::AbstractArray{T},
                                             air_temperature::AbstractArray{T},
                                             wind::AbstractArray{T},
                                             vapour_deficit::AbstractArray{T},
                                             crop_rootdist::AbstractArray{T},
                                             crop_rootzone_available_water::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             CFT::CFTParameters,
                                             kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers, use_precomputed_conductance = kernel_params
    @unpack ALPHAM, GM, LAMBDA_OPT = lpjmlparams
    canopy_height = T(kernel_params.canopy_height)
    # Penman-Monteith's stomatal sensitivity, expressed as LPJmL's own two
    # constants so that no equation downstream changes. `aerodynamic_coupling`
    # is 0 by default, which returns those constants untouched.
    alpha_c, shape_c = coupled_demand_parameters(
        T(ALPHAM), T(GM), pet_eeq[cell], air_temperature[cell],
        vapour_deficit[cell], wind[cell], canopy_height * crop_fpar[cell],
        T(lpjmlparams.aerodynamic_coupling),
    )
    @unpack fpc, emax, gmin, depletion_fraction, depletion_demand_slope = CFT

    co2_index = length(co2) == 1 ? 1 : cell
    if !use_precomputed_conductance
        crop_gp[cell] = compute_canopy_conductance(
            photos_adtmm[cell], co2[co2_index], daylength[cell], crop_fpar[cell],
            T(gmin), T(LAMBDA_OPT),
        )
    end

    wr = zero(T)
    rootzone_water = zero(T)
    for l in 1:soil_layers
        wr += soil_w[l, cell] * crop_rootdist[l]
        if l <= 3
            rootzone_water += soil_w[l, cell] * soil_whcs[l, cell] * crop_rootdist[l]
        end
    end
    crop_rootzone_available_water[cell] = rootzone_water

    if crop_isgrowing[cell] == 1
        # The PT potential, not `demand`: see `demand_adjusted_depletion`. Demand
        # is left computed here rather than below only because the two now read
        # the same `pet_eeq` and keeping them adjacent makes the difference
        # between them visible.
        demand = compute_transpiration_demand(
            crop_canopy_wet[cell], pet_eeq[cell], alpha_c, shape_c, crop_gp[cell],
        )
        adjusted_depletion = demand_adjusted_depletion(
            T(depletion_fraction), pet_eeq[cell] * alpha_c, T(depletion_demand_slope),
        )
        supply = compute_transpiration_supply(
            T(emax), wr, crop_rootc[cell], adjusted_depletion)

        crop_w_demandsum[cell] += demand
        if supply > demand
            crop_w_supplysum[cell] += demand
        else
            crop_w_supplysum[cell] += supply
        end

        crop_wdf[cell] = compute_water_sufficiency(
            crop_w_supplysum[cell], crop_w_demandsum[cell],
        )

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            # The same supply/demand ratio, so it takes the same plateau: `wscal`
            # drives LAI senescence (`lai_crop.jl:52`) and a senescence using a
            # different water stress from allocation would be incoherent.
            crop_wscal[cell] = (emax * compute_available_fraction(wr, adjusted_depletion)) /
                (pet_eeq[cell] * alpha_c / (one(T) + (shape_c * alpha_c) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        # Potential transpiration constrained by demand/supply and canopy fraction.
        if wr > 0
            transp = min(supply, demand) / wr * fpc
        else
            transp = zero(T)
        end

        transp_cor = zero(T)

        # Apply layer-wise extraction cap so uptake does not exceed layer storage.
        if transp > 0
            for l in 1:soil_layers
                transp_tmp, capped = compute_layer_transpiration(
                    transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
                )
                transp_cor += transp_tmp
                capped && transp_cor < T(1e-5) && (transp_cor = zero(T))
            end
        else
            transp_cor = zero(T)
        end

        if wr > 0
            transp = transp_cor / wr
        else
            transp = zero(T)
        end

        # LPJmL recomputes actual canopy conductance after layer extraction.
        # Store it in `gp`; downstream lambda solving consumes this actual value.
        actual_supply = fpc > zero(T) ? transp_cor / fpc : zero(T)
        crop_gp[cell] = compute_actual_canopy_conductance(
            crop_gp[cell], actual_supply, demand, crop_canopy_wet[cell], pet_eeq[cell],
            alpha_c, shape_c,
        )

        # Distribute corrected transpiration back to layers by root distribution.
        for l in 1:soil_layers
            crop_trans_layer[l, cell], _ = compute_layer_transpiration(
                transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
            )
        end
    else
        crop_gp[cell] = zero(T)
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = zero(T)
        end
        crop_w_demandsum[cell] = zero(T)
        crop_w_supplysum[cell] = zero(T)
        crop_wdf[cell] = zero(T)
        # Neutral stress for an absent stand; is_growing still gates all fluxes.
        crop_wscal[cell] = one(T)
    end
end
