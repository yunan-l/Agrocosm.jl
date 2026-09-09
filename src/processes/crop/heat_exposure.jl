# Flowering-window heat exposure, computed without a sub-daily assimilation loop.
#
# The reproductive sink reads one field, `heat_exposure_hours`, and until now
# that field could only be filled by the sub-daily assimilation kernels, which
# accumulate it inside the loop they already run. That made "the sink" and
# "sub-daily photosynthesis" a single choice. They are not, and the measurement
# that separates them is in `docs/07_ablation_framework.md`:
#
#   - The accumulator reads only daily state - the daily mean, the daily range,
#     the daylength, and, through `organ`, humidity, pressure, wind, shortwave,
#     longwave, albedo, LAI and canopy conductance. It never reads anything the
#     assimilation accumulator computes. The two sums merely share a loop.
#     Reimplementing it outside the kernel from a running simulation's daily
#     state reproduces the kernel's value to within 0.04-3.29% at the five gate
#     cells, the residual being the auxiliary state read after the day rather
#     than at kernel time.
#   - Twenty-four evaluations of `organ_leaf_temperature` are cheap. The cost
#     that makes sub-daily integration expensive is the assimilation it wraps:
#     photosynthesis, the nitrogen limitation, and a bisection for Rubisco
#     capacity.
#
# So this file exists to make the cheap half available on its own. But the
# reason to want it is not cost, and this is the part worth stating plainly,
# because it changes what the model claims.
#
# The sub-daily assimilation kernel under-assimilates. That thins the canopy
# (rice LAI 2.82 against 4.92 on the daily kernel), which raises leaf
# temperature (peak 39.1 C against 35.9 C), which raises exposure (141.7 h
# against 73.2 h), which sterilises more grain. The yield deficit against the
# observational reference and the size of the event response are therefore not
# two results, they are one coupled defect. Driving the sink from a trajectory
# whose canopy is right is the honest configuration, not merely the affordable
# one - and the perturbation-B exposure response survives it, rising at rice
# (+31.9% against +12.0%) because a canopy at 39 C has already saturated the
# sterility logistic and has no headroom left to respond with.
#
# One writer, always. `SimulationConfiguration` rejects this switch together
# with `subdaily_photosynthesis`, so `heat_exposure_hours` has exactly one
# source in any run that can be built. Two sources would be worse than either:
# the field would carry whichever kernel ran last, and no ablation cell using it
# would mean anything.
#
# Whether the integral is taken at leaf or at air temperature is still
# `organ_temperature`'s business, exactly as inside the sub-daily kernel: with
# `organ` present `organ_leaf_temperature` solves the canopy energy balance,
# and without it that same function returns the sub-step air temperature. Air
# loses 50-84% of the exposure hours at these cells, so the leaf is
# load-bearing, but keeping the air variant expressible preserves the ablation
# cell that separates the sink mechanism from the departure that triggers it.

"""
    heat_exposure!(CFT, crop, daylength, temperature, diurnal; organ = nothing)

Fill `heat_exposure_hours` with today's duration above the sterility threshold,
integrated over sub-steps of the day without integrating assimilation.

`diurnal === nothing` is a no-op, so the call site can pass the argument
unconditionally and stay type-stable - the same convention `photosynthesis!`
uses. When it is a `DiurnalForcing` the loop is the one the sub-daily
assimilation kernels run for this quantity, with the same integrand, the same
`interval = daylength / steps` convention and the same accumulation order, so
the two agree given the same daily state. `test/processes/crop/test_heat_exposure.jl`
asserts that agreement rather than trusting it.

Runs once per day, before `reproductive_sink!`, which consumes the field.
"""
heat_exposure!(CFT, crop, daylength, temperature, ::Nothing; organ = nothing) = nothing

function heat_exposure!(CFT::CFTParameters,
                        crop,
                        pet_daylength::AbstractArray{T},
                        temp::AbstractArray{T},
                        diurnal::DiurnalForcing;
                        organ = nothing,
) where {T <: AbstractFloat}
    launch_1D!(
        heat_exposure_kernel!,
        crop_stress_auxiliary(crop).heat_exposure_hours,
        crop_stress_auxiliary(crop).filling_exposure_hours,
        pet_daylength,
        temp,
        diurnal.range,
        CFT,
        diurnal.config,
        organ,
    )
    return nothing
end

@kernel inbounds = true function heat_exposure_kernel!(
    heat_exposure_hours::AbstractVector{T},
    filling_exposure_hours::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
    ::DiurnalConfig{STEPS, SHAPE, CAPACITY},
    organ,
) where {T <: AbstractFloat, STEPS, SHAPE, CAPACITY}
    cell = @index(Global)
    @unpack sterility_temperature, filling_temperature = CFT
    @unpack leaf_dimension, leaf_emissivity, lightextcoeff = CFT

    temperature_cell = temperature[cell]
    daylength_cell = daylength[cell]
    range_cell = max(zero(T), diurnal_range[cell])
    interval = daylength_cell / T(STEPS)

    exposure = zero(T)
    filling = zero(T)
    for index in 1:STEPS
        temperature_substep = diurnal_temperature(
            index, STEPS, temperature_cell, range_cell, daylength_cell, SHAPE,
        )
        radiation_fraction = diurnal_radiation_fraction(index, STEPS, SHAPE, T)
        # With organ temperature off this is compile-time dead and the loop
        # integrates duration at sub-step air temperature.
        shortwave_substep = organ === nothing ? zero(T) : diurnal_shortwave_rate(
            radiation_fraction, STEPS, organ.shortwave[cell], daylength_cell,
        )
        leaf_substep = organ_leaf_temperature(
            organ, temperature_substep, cell, shortwave_substep,
            T(leaf_dimension), T(leaf_emissivity), T(lightextcoeff),
        )
        exposure += interval * smooth_exceedance(
            leaf_substep - T(sterility_temperature), T(STERILITY_SMOOTHING_WIDTH),
        )
        # An independent accumulator at the lower grain-filling threshold. It
        # adds one exceedance per sub-step and touches nothing else, so the
        # sterility integral above is unchanged bitwise.
        filling += interval * smooth_exceedance(
            leaf_substep - T(filling_temperature), T(STERILITY_SMOOTHING_WIDTH),
        )
    end
    # Unconditional, matching the sub-daily kernels: a stale value from an
    # earlier day would be read by the sink as if it were today's.
    dark = daylength_cell <= zero(T)
    heat_exposure_hours[cell] = dark ? zero(T) : exposure
    filling_exposure_hours[cell] = dark ? zero(T) : filling
end

"""
    daily_statistic_exposure!(CFT, crop, daylength, temperature, diurnal_range,
                              enabled)

Fill `heat_exposure_hours` from daily aggregates in closed form: no sub-step
loop, no canopy energy balance, air temperature only.

`enabled === false` is a no-op, so the driver can call it unconditionally.

This is the GGCM analogue of the sink's input and it exists to keep an
objection answerable. "A daily-mean kernel cannot see a mean-preserving change
in diurnal range" is true of this model's daily kernel, but the published models
the paper positions itself against read `tasmax` and `tasmin`, so a perturbation
that widens the range raises their maximum. Whether a criterion built from those
two numbers responds is a measurement, and this is the cell that makes it.

It is deliberately the weakest of the three exposure paths - hard threshold,
sinusoid reconstruction, no leaf - because its job is to be the floor the other
two are compared against, not to be good.
"""
daily_statistic_exposure!(CFT, crop, daylength, temperature, diurnal_range,
                          enabled::Bool) = enabled ?
    _daily_statistic_exposure!(CFT, crop, daylength, temperature, diurnal_range) :
    nothing

function _daily_statistic_exposure!(CFT::CFTParameters,
                                    crop,
                                    pet_daylength::AbstractArray{T},
                                    temp::AbstractArray{T},
                                    diurnal_range::AbstractArray{T},
) where {T <: AbstractFloat}
    launch_1D!(
        daily_statistic_exposure_kernel!,
        crop_stress_auxiliary(crop).heat_exposure_hours,
        crop_stress_auxiliary(crop).filling_exposure_hours,
        pet_daylength,
        temp,
        diurnal_range,
        CFT,
    )
    return nothing
end

@kernel inbounds = true function daily_statistic_exposure_kernel!(
    heat_exposure_hours::AbstractVector{T},
    filling_exposure_hours::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack sterility_temperature, filling_temperature = CFT
    range = max(zero(T), diurnal_range[cell])
    heat_exposure_hours[cell] = daily_statistic_exposure_hours(
        temperature[cell], range, daylength[cell], T(sterility_temperature),
    )
    filling_exposure_hours[cell] = daily_statistic_exposure_hours(
        temperature[cell], range, daylength[cell], T(filling_temperature),
    )
end
