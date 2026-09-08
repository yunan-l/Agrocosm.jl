# Opt-in within-day integration of the assimilation core.
#
# Nothing in `photosynthesis.jl` is edited. These kernels sit beside the daily
# ones and are selected by dispatch on the `diurnal` argument, the same pattern
# the weather-AD extension uses. `diurnal === nothing` forwards to the daily
# method unchanged, so production is bitwise unaffected.
#
# What changes: the single evaluation of the instantaneous assimilation rate at
# the daily-mean state becomes a midpoint-rule sum of that same rate over a
# generated diurnal cycle. Concretely, the daily kernel's
#
#     je  = c1 * apar * cmass * cq / (daylength + 1e-5)
#     agd = compute_co_limited_assimilation(je, jc, theta, daylength)
#
# becomes, per sub-step, the same two lines with `apar` replaced by that
# sub-step's share of the day's PAR (its fraction times `steps`) and `daylength`
# replaced by the sub-step length. `compute_co_limited_assimilation` is linear in
# its last argument, so passing the sub-step length gives the sub-step's
# contribution and no new function is needed.
#
# What deliberately does not change:
#
# - Daily bookkeeping. Phenology, allocation, harvest, soil C/N, the water
#   balance, output and checkpointing all remain daily and untouched.
# - Leaf respiration and net assimilation stay daily, gated by the daily
#   temperature stress exactly as today.
# - `comp_vcmax`. The potential-Vcmax branch is an analytic daily-integrated
#   solution; its `sigma` term comes from integrating the daily light curve, so
#   making it sub-daily means re-deriving it. It is left on the daily path.
#   Consequence, stated plainly: on a hot day the capacity `vcmax` is still set
#   from daily-mean stress, so this step corrects the assimilation integral but
#   not the capacity that integral draws on. The sub-daily temperature still
#   reaches the Rubisco-limited branch through `c2`. Making the capacity
#   sub-daily is a separate change with its own validation.
# - `solve_lambda!`, `limit_vcmax_by_nitrogen!` and `recouple_nitrogen_water!`
#   still operate at the daily mean.
#
# Exact degeneracy onto the daily kernel, verified rather than assumed:
#
# | steps | shape             | diurnal range | identical to daily kernel |
# | ----- | ----------------- | ------------- | ------------------------- |
# | 1     | any               | 0             | yes, exactly              |
# | 1     | :flat             | any           | yes, exactly              |
# | 1     | :daytime_neutral  | any           | yes, exactly              |
# | 1     | :sinusoid         | > 0           | no; solar-noon temperature |
#
# The `:daytime_neutral` row holds because its single sub-step subtracts the
# closed-form sub-step mean, which at `steps == 1` is the solar-noon value
# itself.
#
# Every row above is `steps == 1`. With more sub-steps only `:flat` degenerates,
# and then only up to round-off, because summing `steps` contributions is not
# associative in floating point:
#
# | steps | shape     | diurnal range | identical to daily kernel        |
# | ----- | --------- | ------------- | -------------------------------- |
# | > 1   | :flat     | any           | yes, to a few eps                |
# | > 1   | :sinusoid | 0             | no; light curvature alone        |
#
# That last row is easy to misread as a bug and is not one. A zero diurnal
# range removes the temperature spread, but a non-flat shape still distributes
# the day's PAR unevenly across sub-steps, and `compute_co_limited_assimilation`
# is concave in light. Integrating a concave response over an uneven light
# course yields less than evaluating it once at the mean course. So the scheme
# has two independent Jensen channels: temperature curvature, which the diurnal
# range drives, and light curvature, which the shape drives on its own. Only
# `:flat` switches both off.

"""
    DiurnalForcing(config, range)

Host-side pairing of the integration settings with the per-cell diurnal
temperature range (C). Kept separate from `DiurnalConfig` because the config is
a zero-size singleton carrying its settings in the type, while the range is a
backend array; the kernel receives both, the config for free.

`range` is normally `tasmax - tasmin` from the climate forcing. A zero range
reduces every shape to the daily state.
"""
struct DiurnalForcing{C, A}
    config::C
    range::A
end

DiurnalForcing(range; kwargs...) = DiurnalForcing(DiurnalConfig(; kwargs...), range)

diurnal_steps(forcing::DiurnalForcing) = diurnal_steps(forcing.config)
diurnal_shape(forcing::DiurnalForcing) = diurnal_shape(forcing.config)

# `nothing` keeps the existing daily behaviour, so every call site can pass the
# argument unconditionally and stay type-stable.
# `organ` is swallowed here rather than forwarded: leaf temperature is solved
# per sub-step, so it has nowhere to live without the sub-daily loop. The
# configuration layer rejects "organ temperature on, sub-daily off" outright,
# which is where that combination should fail.
photosynthesis!(pathway, CFT, crop, apar, daylength, temperature, co2,
                ::Nothing; organ = nothing, kwargs...) =
    photosynthesis!(pathway, CFT, crop, apar, daylength, temperature, co2; kwargs...)

photosynthesis!(::Val{:C3}, CFT, crop, apar, daylength, temperature, co2,
                diurnal::DiurnalForcing; kwargs...) =
    photosynthesis_subdaily_C3!(CFT, crop, apar, daylength, temperature, co2,
                                diurnal; kwargs...)
photosynthesis!(::Val{:C4}, CFT, crop, apar, daylength, temperature, co2,
                diurnal::DiurnalForcing; kwargs...) =
    photosynthesis_subdaily_C4!(CFT, crop, apar, daylength, temperature,
                                diurnal; kwargs...)

"""Cell-local C3 assimilation with the day's total split over sub-steps."""
function photosynthesis_subdaily_C3!(CFT::CFTParameters,
                                     crop,
                                     apar::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     co2::AbstractArray{T},
                                     diurnal::DiurnalForcing;
                                     lpjmlparams::LPJmLParams = lpjmlparams,
                                     photoparams::PhotoParams = photoparams,
                                     comp_vcmax = false,
                                     organ = nothing,
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_subdaily_c3_kernel!,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_fluxes(crop).carbon.net_assimilation,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_fluxes(crop).carbon.leaf_respiration,
        crop_photosynthesis_auxiliary(crop).potential_vcmax,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        apar,
        pet_daylength,
        temp,
        co2,
        diurnal.range,
        CFT,
        lpjmlparams,
        photoparams,
        comp_vcmax,
        diurnal.config,
        organ,
    )
    return nothing
end

@kernel inbounds = true function photosynthesis_subdaily_c3_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    co2::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
    lpjmlparams::LPJmLParams,
    photoparams::PhotoParams,
    comp_vcmax::Bool,
    ::DiurnalConfig{STEPS, SHAPE},
    organ,
) where {T <: AbstractFloat, STEPS, SHAPE}
    cell = @index(Global)
    @unpack b, path, temp_co2, temp_photos = CFT
    @unpack leaf_dimension, leaf_emissivity, lightextcoeff = CFT
    @unpack ko25, kc25, alphac3, theta, LAMBDA_OPT = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photoparams
    @unpack tmc3, tmc4 = photoparams

    # Daily quantities. These drive the capacity solve and the respiration gate
    # and are identical to the daily kernel.
    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    temperature_cell = temperature[cell]
    daylength_cell = daylength[cell]
    co2_cell = co2[length(co2) == 1 ? 1 : cell]
    range_cell = max(zero(T), diurnal_range[cell])

    ko = T(ko25) * T(q10ko)^((temperature_cell - T(25)) * T(0.1))
    kc = T(kc25) * T(q10kc)^((temperature_cell - T(25)) * T(0.1))
    fac = kc * (one(T) + T(po2) / ko)
    tau = T(tau25) * T(q10tau)^((temperature_cell - T(25)) * T(0.1))
    gammastar = T(po2) / (T(2) * tau)

    if comp_vcmax
        lambda[cell] = T(LAMBDA_OPT)
        if inactive || apar[cell] <= zero(T) || daylength_cell <= zero(T)
            vcmax[cell] = zero(T)
        else
            internal_co2 = T(lambdamc3) * co2_cell
            c1 = stress * T(alphac3) *
                ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
            c2 = (internal_co2 - gammastar) / (internal_co2 + fac)
            s = T(24) / daylength_cell * T(b)
            sigma = one(T) - (c2 - s) / (c2 - T(theta) * s)
            sigma = sqrt(max(zero(T), sigma))
            potential = (one(T) / T(b)) * (c1 / c2) *
                ((T(2) * T(theta) - one(T)) * s -
                 (T(2) * T(theta) * s - c2) * sigma) *
                apar[cell] * T(cmass) * T(cq)
            vcmax[cell] = max(zero(T), potential)
        end
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    # Within-day integral of the instantaneous rate.
    internal_co2 = lambda[cell] * co2_cell
    rubisco_capacity = hour2day(vcmax[cell])
    interval = daylength_cell / T(STEPS)
    gross = zero(T)
    for index in 1:STEPS
        temperature_substep = diurnal_temperature(
            index, STEPS, temperature_cell, range_cell, daylength_cell, SHAPE,
        )
        radiation_fraction = diurnal_radiation_fraction(index, STEPS, SHAPE, T)
        # With organ temperature off this is compile-time dead and the loop is
        # bitwise what step 1 produced.
        shortwave_substep = organ === nothing ? zero(T) : diurnal_shortwave_rate(
            radiation_fraction, STEPS, organ.shortwave[cell], daylength_cell,
        )
        # Every temperature-dependent term below is a leaf process, so all of
        # them follow the leaf, not the air: the stress response and the
        # Michaelis-Menten and specificity constants alike are enzyme kinetics
        # happening inside the leaf.
        leaf_substep = organ_leaf_temperature(
            organ, temperature_substep, cell, shortwave_substep,
            T(leaf_dimension), T(leaf_emissivity), T(lightextcoeff),
        )
        stress_substep = compute_photosynthesis_temperature_stress(
            daylength_cell, leaf_substep, path, temp_co2, temp_photos,
            T(tmc3), T(tmc4),
        )
        ko_substep = T(ko25) * T(q10ko)^((leaf_substep - T(25)) * T(0.1))
        kc_substep = T(kc25) * T(q10kc)^((leaf_substep - T(25)) * T(0.1))
        fac_substep = kc_substep * (one(T) + T(po2) / ko_substep)
        tau_substep = T(tau25) * T(q10tau)^((leaf_substep - T(25)) * T(0.1))
        gammastar_substep = T(po2) / (T(2) * tau_substep)

        c1 = stress_substep * T(alphac3) *
            ((internal_co2 - gammastar_substep) /
             (internal_co2 + T(2) * gammastar_substep))
        c2 = (internal_co2 - gammastar_substep) / (internal_co2 + fac_substep)
        apar_substep = radiation_fraction * T(STEPS) * apar[cell]
        je = c1 * apar_substep * T(cmass) * T(cq) / (daylength_cell + T(1e-5))
        jc = c2 * rubisco_capacity
        agd = compute_co_limited_assimilation(je, jc, T(theta), interval)
        # Gross assimilation cannot be negative within a sub-step, so the floor
        # is applied per sub-step rather than to the daily total. At steps == 1
        # the two are the same expression.
        gross += (stress_substep < T(1e-2)) ? zero(T) : max(zero(T), agd)
    end
    gross_assimilation[cell] = gross

    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    net_assimilation[cell], adt = compute_net_assimilation(gross, leaf, daylength_cell)
    water_limited_assimilation[cell] = compute_water_limited_assimilation(
        adt, T(cmass), temperature_cell, T(p),
    )
end

"""Cell-local C4 assimilation with the day's total split over sub-steps."""
function photosynthesis_subdaily_C4!(CFT::CFTParameters,
                                     crop,
                                     apar::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     diurnal::DiurnalForcing;
                                     lpjmlparams::LPJmLParams = lpjmlparams,
                                     photoparams::PhotoParams = photoparams,
                                     comp_vcmax = false,
                                     organ = nothing,
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_subdaily_c4_kernel!,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_fluxes(crop).carbon.net_assimilation,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_fluxes(crop).carbon.leaf_respiration,
        crop_photosynthesis_auxiliary(crop).potential_vcmax,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        apar,
        pet_daylength,
        temp,
        diurnal.range,
        CFT,
        lpjmlparams,
        photoparams,
        comp_vcmax,
        diurnal.config,
        organ,
    )
    return nothing
end

@kernel inbounds = true function photosynthesis_subdaily_c4_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
    lpjmlparams::LPJmLParams,
    photoparams::PhotoParams,
    comp_vcmax::Bool,
    ::DiurnalConfig{STEPS, SHAPE},
    organ,
) where {T <: AbstractFloat, STEPS, SHAPE}
    cell = @index(Global)
    @unpack b, path, temp_co2, temp_photos = CFT
    @unpack leaf_dimension, leaf_emissivity, lightextcoeff = CFT
    @unpack alphac4, theta, LAMBDA_OPT = lpjmlparams
    @unpack lambdamc4, cmass, cq, p, tmc3, tmc4 = photoparams

    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    temperature_cell = temperature[cell]
    daylength_cell = daylength[cell]
    range_cell = max(zero(T), diurnal_range[cell])

    if comp_vcmax
        lambda[cell] = T(LAMBDA_OPT)
        if inactive || apar[cell] <= zero(T) || daylength_cell <= zero(T)
            vcmax[cell] = zero(T)
        else
            c1 = stress * T(alphac4)
            s = T(24) / daylength_cell * T(b)
            sigma = one(T) - (one(T) - s) / (one(T) - T(theta) * s)
            sigma = sqrt(max(zero(T), sigma))
            potential = (one(T) / T(b)) * c1 *
                ((T(2) * T(theta) - one(T)) * s -
                 (T(2) * T(theta) * s - one(T)) * sigma) *
                apar[cell] * T(cmass) * T(cq)
            vcmax[cell] = max(zero(T), potential)
        end
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    # C4 temperature dependence enters only through the stress scalar, so the
    # sub-step loop is cheaper here than for C3.
    phipi = min(one(T), lambda[cell] / T(lambdamc4))
    rubisco_capacity = hour2day(vcmax[cell])
    interval = daylength_cell / T(STEPS)
    gross = zero(T)
    for index in 1:STEPS
        temperature_substep = diurnal_temperature(
            index, STEPS, temperature_cell, range_cell, daylength_cell, SHAPE,
        )
        radiation_fraction = diurnal_radiation_fraction(index, STEPS, SHAPE, T)
        shortwave_substep = organ === nothing ? zero(T) : diurnal_shortwave_rate(
            radiation_fraction, STEPS, organ.shortwave[cell], daylength_cell,
        )
        leaf_substep = organ_leaf_temperature(
            organ, temperature_substep, cell, shortwave_substep,
            T(leaf_dimension), T(leaf_emissivity), T(lightextcoeff),
        )
        stress_substep = compute_photosynthesis_temperature_stress(
            daylength_cell, leaf_substep, path, temp_co2, temp_photos,
            T(tmc3), T(tmc4),
        )
        c1 = stress_substep * phipi * T(alphac4)
        apar_substep = radiation_fraction * T(STEPS) * apar[cell]
        je = c1 * apar_substep * T(cmass) * T(cq) / (daylength_cell + T(1e-5))
        agd = compute_co_limited_assimilation(je, rubisco_capacity, T(theta), interval)
        gross += (stress_substep < T(1e-2)) ? zero(T) : max(zero(T), agd)
    end
    gross_assimilation[cell] = gross

    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    net_assimilation[cell], adt = compute_net_assimilation(gross, leaf, daylength_cell)
    water_limited_assimilation[cell] = compute_water_limited_assimilation(
        adt, T(cmass), temperature_cell, T(p),
    )
end
