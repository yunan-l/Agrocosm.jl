# The P-model: optimality-based light use efficiency (Prentice et al. 2014,
# Wang et al. 2017, Stocker et al. 2020 GMD 13:1545).
#
# WHY A SECOND PHOTOSYNTHESIS PATH. Measured against a census, this model's yield
# skill IS its assimilation skill: season GPP correlates 0.518 with US county
# maize anomalies and yield 0.538, 0.548 and 0.558 for soybean, 0.260 and 0.280
# for European wheat. Everything between biomass and grain is worth about 4% of
# the correlation (`docs/34`). So assimilation is where 96% of the skill is made
# and the only place a response-function change can move much - and a corrected
# response function is one of only two kinds of change that has ever worked here.
#
# WHAT IS DIFFERENT. LPJmL solves Farquhar with a prescribed lambda and a
# nitrogen-limited Vcmax. The P-model PREDICTS light use efficiency from two
# optimality hypotheses instead of prescribing it:
#
#   * least-cost: the ci/ca ratio minimises the summed cost of maintaining
#     carboxylation and transpiration capacity, which makes it respond to vapour
#     pressure deficit, temperature and pressure without a stomatal model;
#   * coordination: Vcmax acclimates so the Rubisco- and light-limited rates
#     coincide, which removes Vcmax as a free field.
#
# It is therefore an acclimated, weeks-to-months formulation, not an instantaneous
# one - which is the right timescale for a daily crop model and the wrong one for
# an hourly canopy model.
#
# This file is pure functions with no state, so it can be validated against the
# published reference implementation before any of it is wired into the model.

"""Water viscosity relative to its value at 25 C, Huber et al. (2009).

The least-cost ratio carries `1.6 * eta_star` because transpiration cost scales
with viscosity, so this is not a detail: it is what makes the optimal ci/ca
respond to temperature at all.
"""
@inline function pmodel_viscosity_ratio(temperature::T, pressure::T) where {T <: AbstractFloat}
    tk = temperature + T(273.15)
    # Huber's dimensionless density of liquid water, via the Tumlirz equation as
    # `rpmodel` uses it.
    tc = temperature
    lambda = T(1788.316) + T(21.55053) * tc - T(0.4695911) * tc * tc +
             T(3.096363e-3) * tc^3 - T(7.341182e-6) * tc^4
    po = T(5918.499) + T(58.05267) * tc - T(1.1253317) * tc * tc +
         T(6.6123869e-3) * tc^3 - T(1.4661625e-5) * tc^4
    vinf = T(0.6980547) - T(7.435626e-4) * tc + T(3.704258e-5) * tc * tc -
           T(6.315724e-7) * tc^3 + T(9.829576e-9) * tc^4 -
           T(1.197269e-10) * tc^5 + T(1.005461e-12) * tc^6 -
           T(5.437898e-15) * tc^7 + T(1.69946e-17) * tc^8 - T(2.295063e-20) * tc^9
    pbar = pressure * T(1e-5)
    vau = vinf + lambda / (po + pbar)
    rho = T(1000) / vau
    # Huber et al. (2009) viscosity from density and temperature.
    tbar = tk / T(647.096)
    rbar = rho / T(322.0)
    mu0 = T(1e2) * sqrt(tbar) / (T(1.67752) + T(2.20462) / tbar +
          T(0.6366564) / (tbar * tbar) - T(0.241605) / (tbar^3))
    h = (
        (T(0.520094), T(0.0850895), T(-1.08374), T(-0.289555), T(0.0), T(0.0)),
        (T(0.222531), T(0.999115), T(1.88797), T(1.26613), T(0.0), T(0.120573)),
        (T(-0.281378), T(-0.906851), T(-0.772479), T(-0.489837), T(-0.25704), T(0.0)),
        (T(0.161913), T(0.257399), T(0.0), T(0.0), T(0.0), T(0.0)),
        (T(-0.0325372), T(0.0), T(0.0), T(0.0698452), T(0.0), T(0.0)),
        (T(0.0), T(0.0), T(0.0), T(0.0), T(0.00872102), T(0.0)),
        (T(0.0), T(0.0), T(0.0), T(-0.00435673), T(0.0), T(-0.000593264)),
    )
    # Huber's double sum is `sum_i (1/Tbar - 1)^i * sum_j H_ij (rhobar - 1)^j`,
    # and the table above is laid out DENSITY-major - seven rows of six - so the
    # density power indexes the row and the temperature power the column.
    # Swapping the two gives a viscosity of order 1e7 instead of 1e-3, which is
    # how this was caught: against reference water viscosity, not against the
    # P-model output, where it merely shifted xi by 4%.
    ctbar = T(1) / tbar - T(1)
    total = zero(T)
    for i in 1:6
        inner = zero(T)
        for j in 1:7
            inner += h[j][i] * (rbar - one(T))^(j - 1)
        end
        total += ctbar^(i - 1) * inner
    end
    mu1 = exp(rbar * total)
    return mu0 * mu1 * T(1e-6)
end

"""Michaelis-Menten coefficient of Rubisco for carboxylation, Pa.

`Kmm = Kc * (1 + pO2 / Ko)`, both Arrhenius in temperature; oxygen is 20.9476% of
the ambient pressure. Bernacchi et al. (2001).
"""
@inline function pmodel_kmm(temperature::T, pressure::T) where {T <: AbstractFloat}
    tk = temperature + T(273.15)
    arrhenius(value25, activation) =
        value25 * exp(activation * (tk - T(298.15)) / (T(298.15) * T(8.3145) * tk))
    kc = arrhenius(T(39.97), T(79430.0))
    ko = arrhenius(T(27480.0), T(36380.0))
    return kc * (one(T) + T(0.209476) * pressure / ko)
end

"""Photorespiratory compensation point, Pa. Scales with pressure and Arrhenius in
temperature from its 25 C, sea-level value of 4.332 Pa (Bernacchi et al. 2001)."""
@inline function pmodel_gammastar(temperature::T, pressure::T) where {T <: AbstractFloat}
    tk = temperature + T(273.15)
    return T(4.332) * (pressure / T(101325.0)) *
           exp(T(37830.0) * (tk - T(298.15)) / (T(298.15) * T(8.3145) * tk))
end

"""Optimal `ci/ca` from the least-cost hypothesis, with the `xi` that sets it.

    xi  = sqrt(beta * (Kmm + gammastar) / (1.6 * eta_star))
    chi = gammastar/ca + (1 - gammastar/ca) * xi / (xi + sqrt(vpd))

`beta` is the unit cost ratio, 146.0 for C3. The square root in `vpd` is the
signature of the theory: doubling the deficit does not halve the ratio.
"""
@inline function pmodel_optimal_chi(
    temperature::T, pressure::T, vpd::T, ca::T, beta::T,
) where {T <: AbstractFloat}
    kmm = pmodel_kmm(temperature, pressure)
    gammastar = pmodel_gammastar(temperature, pressure)
    eta_star = pmodel_viscosity_ratio(temperature, pressure) /
               pmodel_viscosity_ratio(T(25), pressure)
    deficit = max(vpd, zero(T))
    xi = sqrt(beta * (kmm + gammastar) / (T(1.6) * eta_star))
    ca > zero(T) || return (chi = zero(T), xi = xi, gammastar = gammastar, kmm = kmm)
    chi = gammastar / ca + (one(T) - gammastar / ca) * xi / (xi + sqrt(deficit))
    return (chi = chi, xi = xi, gammastar = gammastar, kmm = kmm)
end

"""The CO2 limitation factor `m`, and `m'` after Jmax limitation.

    m  = (ci - gammastar) / (ci + 2 gammastar)
    m' = m * sqrt(1 - (cstar/m)^(2/3))

`cstar = 0.41`. Below `m = cstar` the root is imaginary and the crop cannot
maintain positive assimilation at that cost ratio, so `m'` is zero rather than
NaN - a guard, not a fudge, and the reason this returns before the power.
"""
@inline function pmodel_limitation(
    chi::T, ca::T, gammastar::T, cstar::T,
) where {T <: AbstractFloat}
    ci = chi * ca
    denominator = ci + T(2) * gammastar
    denominator > zero(T) || return zero(T)
    m = (ci - gammastar) / denominator
    m > cstar || return zero(T)
    return m * sqrt(one(T) - (cstar / m)^(T(2) / T(3)))
end

"""Temperature dependence of the intrinsic quantum yield (Bernacchi et al. 2003),
normalised to its value at 25 C so `kphio` stays the calibrated constant."""
@inline function pmodel_quantum_yield_factor(temperature::T) where {T <: AbstractFloat}
    phi(t) = T(0.352) + T(0.022) * t - T(3.4e-4) * t * t
    reference = phi(T(25))
    reference > zero(T) || return zero(T)
    return max(zero(T), phi(temperature) / reference)
end

"""
    pmodel_light_use_efficiency(temperature, pressure, vpd, ca, kphio, beta, cstar)

Grams of carbon fixed per mole of absorbed photosynthetically active radiation.

`GPP = LUE * absorbed_par`, which is the whole point: light use efficiency is
PREDICTED from temperature, vapour pressure deficit, pressure and CO2 rather than
prescribed. C4 photosynthesis concentrates CO2 at Rubisco and does not
photorespire appreciably, so `m'` is one there and the deficit response vanishes -
pass `cstar = 0` and it falls out of the same expression.
"""
@inline function pmodel_light_use_efficiency(
    temperature::T, pressure::T, vpd::T, ca::T,
    kphio::T, beta::T, cstar::T,
) where {T <: AbstractFloat}
    optimum = pmodel_optimal_chi(temperature, pressure, vpd, ca, beta)
    mprime = cstar > zero(T) ?
             pmodel_limitation(optimum.chi, ca, optimum.gammastar, cstar) : one(T)
    return kphio * pmodel_quantum_yield_factor(temperature) * mprime * T(12.0107)
end

"""Vapour pressure deficit in Pa from specific humidity, pressure and temperature.

The P-model reads a deficit and the forcing carries humidity, so the conversion
belongs beside the model that needs it rather than in a caller.
"""
@inline function pmodel_vapour_pressure_deficit(
    temperature::T, pressure::T, specific_humidity::T,
) where {T <: AbstractFloat}
    saturated = T(611.2) * exp(T(17.62) * temperature / (T(243.12) + temperature))
    q = max(specific_humidity, zero(T))
    actual = q * pressure / (T(0.622) + T(0.378) * q)
    return max(saturated - actual, zero(T))
end
