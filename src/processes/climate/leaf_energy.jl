# Canopy energy balance and the leaf-air temperature departure.
#
# Every function here is scalar, allocation-free and branch-light so it is safe
# to call from a CPU or GPU kernel, and smooth in its floating-point arguments
# so the Enzyme reverse pass stays finite. Nothing here iterates: the departure
# is a closed form evaluated twice, which is what keeps the reverse pass
# straight-line code. See docs/04_organ_temperature_design.md.
#
# The physical content is one equation. Balancing absorbed radiation against
# sensible and latent loss, and linearising both the outgoing longwave and the
# saturation vapour pressure about an expansion point, gives
#
#   dT = [Rn - (lambda/p) gv D] / [cp gH + 4 eps sigma T^3 + (lambda/p) gv s]
#
# Every term in the denominator is positive, so it never vanishes and the
# solution is bounded for any input.
#
# Two properties are worth stating because they are what the tests check:
#
#   - As gH grows the departure vanishes: a leaf glued to the air by turbulence
#     is at air temperature. This is the exact degeneracy that makes the
#     ablation experiment's "organ temperature off" arm exact rather than
#     nominal.
#   - The departure changes sign. With stomata open and D modest, the latent
#     term dominates and the canopy sits BELOW air temperature - transpirational
#     cooling. As conductance falls, that term collapses and the canopy runs
#     hot. Reproducing both signs from one equation is the point: measured
#     irrigated rye sits near -2 C while rainfed rye on sandy soil reaches
#     +7.5 C (Siebert et al. 2014, ERL 9:044012).
#
# Saturation vapour pressure deliberately reuses the coefficients already
# embedded in `compute_equilibrium_evaporation` in radiation.jl, where the
# constant 2.503e6 is exactly 610.78 * 17.269 * 237.3. A model must not carry
# two different saturation curves.
#
# UNITS. Every temperature crossing this module's interface is in **degrees
# Celsius**, matching `dailyWeather.temp` and the -100..70 range
# `_validate_climate(:temp, ...)` enforces. Absolute temperature is formed
# locally, and only where the physics demands it: the Stefan-Boltzmann term and
# the ideal-gas density. The returned departure is a temperature *difference*,
# so degrees Celsius and kelvin coincide and no conversion is owed to the
# caller. The remaining units are SI throughout - Pa for pressures, W m-2 for
# fluxes, mol m-2 s-1 for conductances - with the single exception of
# `stomatal_conductance`, which arrives in mm s-1 because that is how
# `crop_canopy_auxiliary(crop).canopy_conductance` is stored, and is converted
# on entry.

const STEFAN_BOLTZMANN = 5.670374419e-8      # W m-2 K-4
const MOLAR_HEAT_CAPACITY_AIR = 29.3         # J mol-1 K-1
const MOLAR_GAS_CONSTANT = 8.314462618       # J mol-1 K-1
const WATER_AIR_MOLAR_RATIO = 0.622          # M_w / M_air
const KELVIN_OFFSET = 273.15
const SATURATION_REFERENCE = 610.78          # Pa, e_s at 0 C
const SATURATION_SLOPE_COEFFICIENT = 17.269  # Tetens, as used by LPJmL
const SATURATION_TEMPERATURE_OFFSET = 237.3  # C, as used by LPJmL

"""
    saturation_vapour_pressure(temperature)

Saturation vapour pressure (Pa) at `temperature` (C), on the same Tetens curve
LPJmL's equilibrium evaporation already uses.
"""
@inline function saturation_vapour_pressure(temperature::T) where {T <: AbstractFloat}
    offset = T(SATURATION_TEMPERATURE_OFFSET) + temperature
    return T(SATURATION_REFERENCE) *
        exp(T(SATURATION_SLOPE_COEFFICIENT) * temperature / offset)
end

"""
    saturation_vapour_pressure_slope(temperature)

`d(e_s)/dT` (Pa K-1) at `temperature` (C), the analytic derivative of
[`saturation_vapour_pressure`](@ref). Identical to the `slope` term inside
`compute_equilibrium_evaporation`.
"""
@inline function saturation_vapour_pressure_slope(temperature::T) where {T <: AbstractFloat}
    offset = T(SATURATION_TEMPERATURE_OFFSET) + temperature
    return T(SATURATION_REFERENCE * SATURATION_SLOPE_COEFFICIENT *
             SATURATION_TEMPERATURE_OFFSET) *
        exp(T(SATURATION_SLOPE_COEFFICIENT) * temperature / offset) / (offset * offset)
end

"""
    vapour_pressure_from_specific_humidity(specific_humidity, pressure)

Actual vapour pressure (Pa) from specific humidity (kg water per kg moist air)
and surface pressure (Pa).

Specific humidity rather than relative humidity is the right sub-daily input:
it is near conserved within a day, so it stays constant across sub-steps while
`e_s(T)` tracks the sub-step temperature, and the vapour pressure deficit peaks
at solar noon on its own.
"""
@inline function vapour_pressure_from_specific_humidity(
    specific_humidity::T, pressure::T,
) where {T <: AbstractFloat}
    q = max(zero(T), specific_humidity)
    return q * pressure / (T(WATER_AIR_MOLAR_RATIO) + (one(T) - T(WATER_AIR_MOLAR_RATIO)) * q)
end

"""
    molar_air_density(temperature, pressure)

Molar density of air (mol m-3) at `temperature` (C) and `pressure` (Pa), used to
convert the model's `mm s-1` conductances onto the `mol m-2 s-1` basis the
energy balance works in.
"""
@inline function molar_air_density(temperature::T, pressure::T) where {T <: AbstractFloat}
    return pressure / (T(MOLAR_GAS_CONSTANT) * (temperature + T(KELVIN_OFFSET)))
end

"""
    latent_heat_of_vaporisation(temperature)

Molar latent heat of vaporisation (J mol-1) at `temperature` (C). Converted from
the mass-basis expression already used in `compute_equilibrium_evaporation`
(`2.495e6 - 2380 T`, J kg-1) so the two stay consistent.
"""
@inline function latent_heat_of_vaporisation(temperature::T) where {T <: AbstractFloat}
    return (T(2.495e6) - T(2380) * temperature) * T(0.018015)
end

"""
    leaf_boundary_layer_conductance(wind_speed, leaf_dimension)

One-sided leaf boundary-layer conductance to heat (mol m-2 s-1) under forced
convection, with the customary 1.4 outdoor-turbulence enhancement:
`1.4 * 0.135 * sqrt(u / d)`.

The 0.135 coefficient is already on a molar basis - it follows from
`Nu = 0.664 Re^0.5 Pr^(1/3)` together with the molar density of air at standard
conditions - so the result must **not** be multiplied by `molar_air_density`
again. Its residual temperature and pressure dependence is well under the
scheme's other biases and is ignored.

`leaf_dimension` is the characteristic dimension in the direction of flow (m),
not the leaf width. Free convection is deliberately omitted: the crossover wind
speed is 0.03-0.17 m s-1 depending on leaf size, and above 1 m s-1 a mixed
formulation differs by under 1%. `wind_speed` is floored rather than allowed to
reach zero, which both keeps the square root differentiable and stands in for
the free-convection floor.
"""
@inline function leaf_boundary_layer_conductance(
    wind_speed::T, leaf_dimension::T,
) where {T <: AbstractFloat}
    speed = max(wind_speed, T(0.1))
    dimension = max(leaf_dimension, T(1e-3))
    return T(1.4) * T(0.135) * sqrt(speed / dimension)
end

"""
    effective_leaf_area(lai, extinction)

Beer-Lambert weighted leaf area `(1 - exp(-k L)) / k` (m2 m-2) that actually
participates in turbulent exchange. Wind and turbulence decay through the
canopy, so the full `lai` overestimates the exchanging area.
"""
@inline function effective_leaf_area(lai::T, extinction::T) where {T <: AbstractFloat}
    k = max(extinction, T(1e-3))
    return (one(T) - exp(-k * max(lai, zero(T)))) / k
end

"""
    canopy_cover_fraction(lai, extinction)

Fraction of incident radiation intercepted by the canopy, `1 - exp(-k L)`. A
sparse canopy intercepts less energy and therefore departs less from air
temperature; at `lai == 0` this is zero and the departure vanishes exactly.
"""
@inline function canopy_cover_fraction(lai::T, extinction::T) where {T <: AbstractFloat}
    return one(T) - exp(-max(extinction, T(1e-3)) * max(lai, zero(T)))
end

"""
    series_conductance(stomatal, boundary)

Two conductances in series, `1/(1/a + 1/b)`, written as `a*b/(a+b)` so a zero
conductance gives zero rather than a division by zero.
"""
@inline function series_conductance(stomatal::T, boundary::T) where {T <: AbstractFloat}
    total = stomatal + boundary
    return total <= zero(T) ? zero(T) : stomatal * boundary / total
end

"""
    _departure_about(expansion_temperature, air_temperature, ...)

Solve the linearised energy balance for `T_leaf - T_air`, expanding the
saturation curve and the outgoing longwave about `expansion_temperature` rather
than about the air temperature.

With `x = T_leaf - T_air` and `m = T_expansion - T_air`, the balance reduces to

    x = [Rn* - (lambda/p) gv D* + m (4 eps sigma T*^3 + (lambda/p) gv s*)]
        / [cp gH + 4 eps sigma T*^3 + (lambda/p) gv s*]

which at `m == 0` is the ordinary first-order form. Calling this twice, the
second time with the first estimate as the expansion point, is algebraically a
single Newton correction and cuts the worst-case linearisation error from
0.46 K to 0.003 K.
"""
@inline function _departure_about(
    expansion_temperature::T,
    air_temperature::T,
    shortwave_absorbed::T,
    net_longwave::T,
    vapour_pressure_air::T,
    pressure::T,
    vapour_conductance::T,
    heat_conductance::T,
    emissivity::T,
    cover_fraction::T,
) where {T <: AbstractFloat}
    offset = expansion_temperature - air_temperature
    expansion_kelvin = expansion_temperature + T(KELVIN_OFFSET)
    air_kelvin = air_temperature + T(KELVIN_OFFSET)
    longwave_coefficient = T(4) * emissivity * T(STEFAN_BOLTZMANN) *
        expansion_kelvin * expansion_kelvin * expansion_kelvin

    # `net_longwave` is the isothermal net longwave the climate forcing supplies
    # (negative for a loss), evaluated at air temperature. Moving the expansion
    # point off the air temperature changes the reference emission, and that
    # shift has to be taken out here in full non-linear form; only the residual
    # departure from the expansion point is left to the linear
    # `longwave_coefficient` term in the denominator. On the first pass
    # `expansion_temperature == air_temperature` and this correction is exactly
    # zero -- which is why omitting it still gave a plausible first estimate but
    # silently destroyed the accuracy the re-expansion is supposed to buy, since
    # the latent term does track the expansion point through `es(T*)`.
    net_radiation = cover_fraction * (
        shortwave_absorbed + net_longwave -
        emissivity * T(STEFAN_BOLTZMANN) *
        (expansion_kelvin^4 - air_kelvin^4)
    )

    latent_coefficient = latent_heat_of_vaporisation(air_temperature) *
        vapour_conductance / pressure
    deficit = saturation_vapour_pressure(expansion_temperature) - vapour_pressure_air
    slope = saturation_vapour_pressure_slope(expansion_temperature)

    damping = cover_fraction * longwave_coefficient + latent_coefficient * slope
    denominator = T(MOLAR_HEAT_CAPACITY_AIR) * heat_conductance + damping
    numerator = net_radiation - latent_coefficient * deficit + offset * damping
    return numerator / denominator
end

"""
    leaf_temperature_departure(air_temperature, shortwave_absorbed, net_longwave,
                               specific_humidity, pressure, wind_speed,
                               stomatal_conductance, lai, leaf_dimension,
                               emissivity, extinction)

Leaf-minus-air temperature (K) for one sub-step, from the closed-form linearised
energy balance with a single re-expansion.

`stomatal_conductance` is the bulk canopy conductance in `mm s-1` per unit
ground area, as stored in `crop_canopy_auxiliary(crop).canopy_conductance`. The
boundary-layer conductance is per unit leaf area and is raised to canopy scale
by twice the effective leaf area, because a flat leaf exchanges sensible heat
from both faces; omitting that factor is worth up to 50% of the flux
(Schymanski & Or 2017, HESS 21:685-706). Longwave is not scaled by leaf area:
the canopy radiates to the sky as a surface.
"""
@inline function leaf_temperature_departure(
    air_temperature::T,
    shortwave_absorbed::T,
    net_longwave::T,
    specific_humidity::T,
    pressure::T,
    wind_speed::T,
    stomatal_conductance::T,
    lai::T,
    leaf_dimension::T,
    emissivity::T,
    extinction::T,
) where {T <: AbstractFloat}
    cover = canopy_cover_fraction(lai, extinction)
    cover <= zero(T) && return zero(T)

    density = molar_air_density(air_temperature, pressure)
    leaf_area = effective_leaf_area(lai, extinction)
    # Both faces of every leaf exchange sensible heat and water vapour.
    two_sided = T(2) * leaf_area
    boundary_heat = two_sided *
        leaf_boundary_layer_conductance(wind_speed, leaf_dimension)
    # Vapour transfer is slightly more efficient than heat over the same layer.
    boundary_vapour = boundary_heat * T(0.147 / 0.135)

    stomatal_molar = max(stomatal_conductance, zero(T)) * T(1e-3) * density
    vapour_conductance = series_conductance(stomatal_molar, boundary_vapour)
    vapour_pressure_air = vapour_pressure_from_specific_humidity(specific_humidity, pressure)

    first_estimate = _departure_about(
        air_temperature, air_temperature, shortwave_absorbed, net_longwave,
        vapour_pressure_air, pressure, vapour_conductance, boundary_heat,
        emissivity, cover,
    )
    # One re-expansion about the first estimate. Fixed cost, no branch, and it
    # is what makes the linearisation error negligible beside the scheme's
    # structural biases rather than comparable to them.
    return _departure_about(
        air_temperature + first_estimate, air_temperature, shortwave_absorbed,
        net_longwave, vapour_pressure_air, pressure, vapour_conductance,
        boundary_heat, emissivity, cover,
    )
end

"""
    OrganTemperatureForcing(specific_humidity, pressure, wind, net_longwave,
                            albedo, lai, conductance)

Per-cell fields the canopy energy balance needs beyond what the sub-daily
assimilation kernel already carries, bundled so the kernel gains one argument
rather than seven. Same pattern as [`DiurnalForcing`](@ref): all fields are
device arrays indexed by cell, and the whole object is `Enzyme.Const` because
none of it is a differentiated control.

`shortwave` is the model's 24-hour-mean `swr` and must never be used as-is
inside the day: the kernel holds the sub-step radiation weights and converts it
with [`diurnal_shortwave_rate`](@ref) first.

Two type parameters rather than one: the humidity and pressure rows arrive as
`view`s into the `(day, cell)` climate matrices, while the remaining fields are
the weather and canopy arrays themselves. A single parameter would demand they
all share a type and reject the pairing outright.
"""
struct OrganTemperatureForcing{V, A}
    specific_humidity::V  # kg kg-1
    pressure::V           # Pa
    wind::A               # m s-1
    shortwave::A          # W m-2, 24-hour mean; convert before use
    net_longwave::A       # W m-2, isothermal, negative for a loss
    albedo::A             # 0-1
    lai::A                # m2 m-2
    conductance::A        # mm s-1, bulk canopy
end

"""
    organ_leaf_temperature(organ, air_temperature, cell, shortwave, dimension,
                           emissivity, extinction)

Sub-step leaf temperature (C). With `organ === nothing` this returns the air
temperature unchanged, and the branch is resolved at compile time, so the
switched-off path is bitwise identical to step 1 rather than merely close.
"""
@inline organ_leaf_temperature(
    ::Nothing, air_temperature::T, cell::Integer, shortwave::T, dimension::T,
    emissivity::T, extinction::T,
) where {T <: AbstractFloat} = air_temperature

@inline function organ_leaf_temperature(
    organ::OrganTemperatureForcing, air_temperature::T, cell::Integer,
    shortwave::T, dimension::T, emissivity::T, extinction::T,
) where {T <: AbstractFloat}
    absorbed = (one(T) - organ.albedo[cell]) * shortwave
    return air_temperature + leaf_temperature_departure(
        air_temperature, absorbed, organ.net_longwave[cell],
        organ.specific_humidity[cell], organ.pressure[cell], organ.wind[cell],
        organ.conductance[cell], organ.lai[cell], dimension, emissivity, extinction,
    )
end

"""
    diurnal_shortwave_rate(fraction, steps, daily_shortwave, daylength)

Instantaneous shortwave irradiance (W m-2) during a sub-step, from the model's
24-hour-mean `swr` (also W m-2) and that sub-step's share of the day's
radiation.

The day's energy is `swr * 86400`; a sub-step receives `fraction` of it over
`daylength * 3600 / steps` seconds, so the rate is
`swr * 24 * fraction * steps / daylength`. Under `:flat` this reduces to
`swr * 24 / daylength`, the whole day's mean power concentrated into the
daylight window, which is the expected sanity check.
"""
@inline function diurnal_shortwave_rate(
    fraction::T, steps::Integer, daily_shortwave::T, daylength::T,
) where {T <: AbstractFloat}
    return daily_shortwave * T(24) * fraction * T(steps) / (daylength + T(1e-5))
end

"""
    leaf_energy_residual(leaf_temperature, air_temperature, ...)

Residual (W m-2) of the **full non-linear** energy balance at `leaf_temperature`:
absorbed radiation minus sensible minus latent loss, with no linearisation
anywhere. Zero at the exact solution.

This exists for the tests. Gate E3 substitutes the closed-form solution back
into this and requires the residual to stay small across the extreme corner of
the parameter space; gate E4 requires the re-expansion to reduce it by at least
two orders of magnitude.
"""
@inline function leaf_energy_residual(
    leaf_temperature::T,
    air_temperature::T,
    shortwave_absorbed::T,
    net_longwave::T,
    specific_humidity::T,
    pressure::T,
    wind_speed::T,
    stomatal_conductance::T,
    lai::T,
    leaf_dimension::T,
    emissivity::T,
    extinction::T,
) where {T <: AbstractFloat}
    cover = canopy_cover_fraction(lai, extinction)
    density = molar_air_density(air_temperature, pressure)
    leaf_area = effective_leaf_area(lai, extinction)
    two_sided = T(2) * leaf_area
    boundary_heat = two_sided *
        leaf_boundary_layer_conductance(wind_speed, leaf_dimension)
    boundary_vapour = boundary_heat * T(0.147 / 0.135)
    stomatal_molar = max(stomatal_conductance, zero(T)) * T(1e-3) * density
    vapour_conductance = series_conductance(stomatal_molar, boundary_vapour)
    vapour_pressure_air = vapour_pressure_from_specific_humidity(specific_humidity, pressure)

    # `net_longwave` is the isothermal value, at air temperature. The leaf's
    # departure changes its emission by the full non-linear T^4 difference --
    # no linearisation anywhere in this function, which is the whole point of
    # having it. At `leaf_temperature == air_temperature` this reduces to
    # `net_longwave` exactly.
    leaf_kelvin = leaf_temperature + T(KELVIN_OFFSET)
    air_kelvin = air_temperature + T(KELVIN_OFFSET)
    net_radiation = cover * (
        shortwave_absorbed + net_longwave -
        emissivity * T(STEFAN_BOLTZMANN) * (leaf_kelvin^4 - air_kelvin^4)
    )
    sensible = T(MOLAR_HEAT_CAPACITY_AIR) * boundary_heat *
        (leaf_temperature - air_temperature)
    latent = latent_heat_of_vaporisation(air_temperature) * vapour_conductance *
        (saturation_vapour_pressure(leaf_temperature) - vapour_pressure_air) / pressure
    return net_radiation - sensible - latent
end
