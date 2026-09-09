# Within-day (diurnal) profiles and integration weights for the assimilation
# core. Every function here is scalar, allocation-free and branch-light so it is
# safe to call from a CPU or GPU kernel, and smooth in its floating-point
# arguments so the Enzyme reverse pass stays finite.
#
# The discrete weights are normalised, not merely sampled from a continuous
# profile. Midpoint sampling of a half-sine does not conserve its integral at
# finite step count (N = 1 overshoots by pi/2), so the daily totals are imposed
# on the weights themselves. This is what makes the daily-mean and sub-daily
# runs share identical daily forcing and identical daily totals, and it makes
# N = 1 degenerate exactly onto the existing daily kernel.
#
# Shape codes are integers, not Symbols, so they can cross a GPU kernel and stay
# constant at the Enzyme boundary:
#
#   1 = :flat            no diurnal cycle; every sub-step sees the daily state
#   2 = :sinusoid        mean-conserving over the full 24 h (physical default)
#   3 = :daytime_neutral sinusoid shifted so the sub-step mean is the daily mean
#
# Shape 3 exists to separate the two things that switching on a diurnal cycle
# does at once: it exposes the crop to curvature within the day (Jensen), and it
# raises the mean temperature seen during daylight above the daily mean. Shape 3
# keeps only the first; shape 2 has both. Reporting them separately is what makes
# the temporal-aggregation result interpretable rather than a single lumped
# number.

const DIURNAL_FLAT = 1
const DIURNAL_SINUSOID = 2
const DIURNAL_DAYTIME_NEUTRAL = 3

"""
    diurnal_substep_time(index, steps, daylength)

Local solar time (h) of the midpoint of sub-step `index` of `steps` equal
sub-intervals spanning the daylight window, which is centred on solar noon.
"""
@inline function diurnal_substep_time(index::Integer, steps::Integer, daylength::T) where {T <: AbstractFloat}
    sunrise = T(12) - daylength * T(0.5)
    return sunrise + (T(index) - T(0.5)) * daylength / T(steps)
end

"""
    diurnal_radiation_weight(index, steps)

Fraction of the day's absorbed PAR received in sub-step `index` of `steps`.
Weights follow a half-sine over the daylight window and are normalised so that
`sum(weights) == 1` analytically for every `steps`, using
`sum_{i=1}^{N} sin(pi (i - 1/2) / N) = 1 / sin(pi / (2N))`.

At `steps == 1` the weight is exactly one, so the sub-step rate becomes
`daily_apar / daylength`: the quantity the existing daily kernel uses.
"""
@inline function diurnal_radiation_weight(index::Integer, steps::Integer, ::Type{T}) where {T <: AbstractFloat}
    steps == 1 && return one(T)
    phase = T(pi) * (T(index) - T(0.5)) / T(steps)
    return sin(phase) * sin(T(pi) / (T(2) * T(steps)))
end

"""
    DiurnalConfig(; steps, shape)

Runtime selection of the within-day integration, carried in the **type** rather
than in fields.

`steps` and `shape` are type parameters, so `DiurnalConfig` is a zero-size
singleton: it crosses a GPU kernel for free, it is trivially constant at the
Enzyme boundary, and inside a kernel the sub-step trip count and the shape
branch are compile-time constants.

That last point is the reason for this design. With a runtime trip count,
reverse-mode AD has to tape dynamic control flow: 48 sub-steps over a 150-day
season is roughly 7,000 taped iterations per cell, which is where Enzyme becomes
slow or its type analysis gives up. A statically known count lets LLVM unroll the
loop and leaves Enzyme straight-line code. Nothing here carries a scientific
coefficient, so it does not belong in `LPJmLParams`.

`capacity_optimum` selects how Rubisco capacity is set when the assimilation
kernel computes it (`comp_vcmax`). The default `false` keeps the Haxeltine-
Prentice analytic solution, whose `sigma` term integrates a *flat* day, exactly
as the daily kernel uses it. `true` re-solves the same optimality problem
against the sub-daily light course the kernel actually integrates
(`subdaily_optimal_vcmax`).

`true` is the internally consistent choice within the optimality hypothesis and
raises capacity by about 1.3x, but it is NOT the default, for three reasons
measured on real cells:

  * The flat-day optimum is not an approximation inside LPJmL, it is the
    definition of `vcmax` that leaf nitrogen demand, `b` and every crop
    calibration coefficient were tuned against. Re-solving it rescales the
    nitrogen cycle away from that calibration, and the four-crop check moved
    wheat's yield the wrong way (-20.9% -> -23.2% sub-daily penalty) while
    leaving soybean untouched.
  * At re-optimised capacity the sub-daily *gross* assimilation is no longer
    bounded above by the flat-day value, so the module's headline mechanism -
    resolving the day lowers assimilation - stops being provable.
  * The solve is a 40-iteration bisection in the innermost kernel and appears
    on the AD tape.

Off it is compile-time dead: the branch disappears and the kernel is bitwise
what the flat-day capacity produces. On, it is the sensitivity experiment that
answers "is your sub-daily loss an artifact of holding capacity at the flat-day
optimum?" with a number - it recovers about 4 pp of rice's 32% penalty and none
of soybean's.

`steps == 1` degenerates onto the daily kernel exactly only for `:flat` and
`:daytime_neutral`, or for any shape at a zero diurnal range. With `:sinusoid`
and a non-zero range it does NOT: the single sub-step samples the midpoint of
the daylight window, i.e. the solar-noon temperature, not the daily mean. The
authoritative table is in `photosynthesis_subdaily.jl`. Measured on a Michigan
soybean cell, `steps = 1, :sinusoid` costs 33% of yield against the daily
kernel, and the step sweep is non-monotone below 8 (-33/-16/-41/-47/-49/-49 at
1/2/3/4/6/8), converging by 8. Treat `steps = 1` as a degeneracy probe, not as
a way to switch the scheme down.
"""
struct DiurnalConfig{STEPS, SHAPE, CAPACITY} end

function DiurnalConfig(; steps::Integer = 1, shape::Integer = DIURNAL_SINUSOID,
                       capacity_optimum::Bool = false)
    steps >= 1 || throw(ArgumentError("steps must be at least 1"))
    shape in (DIURNAL_FLAT, DIURNAL_SINUSOID, DIURNAL_DAYTIME_NEUTRAL) ||
        throw(ArgumentError("shape must be DIURNAL_FLAT, DIURNAL_SINUSOID or DIURNAL_DAYTIME_NEUTRAL"))
    return DiurnalConfig{Int(steps), Int(shape), capacity_optimum}()
end

diurnal_steps(::DiurnalConfig{STEPS}) where {STEPS} = STEPS
diurnal_shape(::DiurnalConfig{STEPS, SHAPE}) where {STEPS, SHAPE} = SHAPE
diurnal_capacity_optimum(
    ::DiurnalConfig{STEPS, SHAPE, CAPACITY},
) where {STEPS, SHAPE, CAPACITY} = CAPACITY

"""
    diurnal_shape_code(shape::Symbol)

Map the configuration symbol onto the integer code. A `Symbol` cannot cross a
GPU kernel, so this runs once on the host.
"""
diurnal_shape_code(shape::Symbol) =
    shape === :flat ? DIURNAL_FLAT :
    shape === :sinusoid ? DIURNAL_SINUSOID :
    shape === :daytime_neutral ? DIURNAL_DAYTIME_NEUTRAL :
    throw(ArgumentError("diurnal shape must be :flat, :sinusoid or :daytime_neutral, got $shape"))

"""
    diurnal_radiation_fraction(index, steps, shape, T)

Fraction of the day's absorbed PAR in sub-step `index`, summing to one over the
sub-steps for every shape and every `steps`.

The assimilation kernel multiplies this by `steps` and then divides by
`daylength + 1e-5`, reusing the existing guard, rather than calling
`diurnal_apar_rate`. That is what makes `steps == 1` reproduce the daily
expression bit for bit instead of only to round-off.
"""
@inline function diurnal_radiation_fraction(
    index::Integer, steps::Integer, shape::Integer, ::Type{T},
) where {T <: AbstractFloat}
    shape == DIURNAL_FLAT && return one(T) / T(steps)
    return diurnal_radiation_weight(index, steps, T)
end

"""
    diurnal_apar_rate(index, steps, daily_apar, daylength, shape)

Instantaneous absorbed PAR (same units as `daily_apar` divided by hours) during
sub-step `index`. For every shape and every `steps`,
`sum_i rate_i * (daylength / steps) == daily_apar` holds by construction, so no
radiation is created or destroyed by refining the sub-daily step.
"""
@inline function diurnal_apar_rate(
    index::Integer, steps::Integer, daily_apar::T, daylength::T, shape::Integer,
) where {T <: AbstractFloat}
    daylength <= zero(T) && return zero(T)
    weight = shape == DIURNAL_FLAT ? one(T) / T(steps) :
        diurnal_radiation_weight(index, steps, T)
    return weight * daily_apar * T(steps) / daylength
end

"""
    diurnal_shape_value(index, steps, daylength)

Unit-amplitude sinusoid `sin(2 pi (t - 8) / 24)` at the sub-step midpoint. It
peaks at 14:00 and troughs at 02:00 local solar time.
"""
@inline function diurnal_shape_value(index::Integer, steps::Integer, daylength::T) where {T <: AbstractFloat}
    time = diurnal_substep_time(index, steps, daylength)
    return sin(T(2) * T(pi) * (time - T(8)) / T(24))
end

"""
    diurnal_shape_substep_mean(steps, daylength)

Mean of `diurnal_shape_value` over the `steps` sub-step midpoints, evaluated in
closed form so it costs no loop inside a kernel.

With `t_i = 12 - L/2 + (i - 1/2) L / N` the sub-step phases form an arithmetic
progression, so the Dirichlet-kernel identity gives

    (1/N) sum_i sin(a + (i - 1/2) d) = sin(a + N d / 2) * sin(N d / 2) /
                                       (N sin(d / 2) * ... )

which reduces to the expression below with `a0` the phase at solar noon and
`half = pi L / 24`. At `steps == 1` it returns the single midpoint value, the
sinusoid at solar noon.
"""
@inline function diurnal_shape_substep_mean(steps::Integer, daylength::T) where {T <: AbstractFloat}
    noon_phase = T(2) * T(pi) * (T(12) - T(8)) / T(24)   # phase at solar noon
    half = T(pi) * daylength / T(24)                      # half-window in phase
    steps == 1 && return sin(noon_phase)
    step_phase = T(2) * half / T(steps)
    denominator = T(steps) * sin(step_phase * T(0.5))
    # The progression is symmetric about solar noon, so the sum collapses to a
    # single sine times a Dirichlet factor. Guard the degenerate zero window.
    abs(denominator) < eps(T) && return sin(noon_phase)
    return sin(noon_phase) * sin(half) / denominator
end

"""
    diurnal_temperature(index, steps, mean_temperature, amplitude, daylength, shape)

Air temperature (C) during sub-step `index`.

- `DIURNAL_FLAT` returns `mean_temperature` unchanged.
- `DIURNAL_SINUSOID` adds `amplitude / 2` times the unit sinusoid. Its mean over
  the full 24 h is `mean_temperature` exactly, for every amplitude.
- `DIURNAL_DAYTIME_NEUTRAL` removes the closed-form sub-step mean of the shape,
  so `(1/N) sum_i T_i == mean_temperature` exactly, for every `steps` and every
  amplitude.

`amplitude` is the diurnal temperature range in C.
"""
@inline function diurnal_temperature(
    index::Integer, steps::Integer, mean_temperature::T, amplitude::T,
    daylength::T, shape::Integer,
) where {T <: AbstractFloat}
    shape == DIURNAL_FLAT && return mean_temperature
    value = diurnal_shape_value(index, steps, daylength)
    shape == DIURNAL_DAYTIME_NEUTRAL &&
        (value -= diurnal_shape_substep_mean(steps, daylength))
    return mean_temperature + amplitude * T(0.5) * value
end

"""
    smooth_exceedance(excess, width)

Logistic replacement for the indicator `excess > 0`, written so neither tail
overflows. `width -> 0` recovers the hard threshold; a hard indicator inside the
day would put a non-differentiable kink in the reverse pass.
"""
@inline function smooth_exceedance(excess::T, width::T) where {T <: AbstractFloat}
    width <= zero(T) && return excess > zero(T) ? one(T) : zero(T)
    scaled = excess / width
    return scaled >= zero(T) ? one(T) / (one(T) + exp(-scaled)) :
        exp(scaled) / (one(T) + exp(scaled))
end

"""
    daily_statistic_exposure_hours(mean_temperature, amplitude, daylength,
                                   critical_temperature)

Daylight hours above `critical_temperature`, in closed form from daily
aggregates alone - no sub-step loop, no energy balance.

This is the GGCM analogue of the sink's input, and it exists to answer an
objection rather than to be the best available physics. The claim that a
daily-mean kernel cannot see a mean-preserving change in diurnal range is true
of *this model's* daily kernel, but it is not true of the published models the
paper compares itself to: they read `tasmax` and `tasmin`, so a perturbation
that raises the range raises their maximum too. A criterion built from those two
numbers therefore has to be measured, not assumed inert.

The temperature course is the same one `diurnal_temperature` reconstructs,

    T(t) = mean + (amplitude / 2) * sin(pi * (t - 8) / 12)

so `T > critical` holds on `t in (8 + 12 asin(z) / pi, 20 - 12 asin(z) / pi)`
with `z = (critical - mean) / (amplitude / 2)`, intersected with the daylight
window `[12 - daylength/2, 12 + daylength/2]`. That is exact, not a quadrature,
and it costs one `asin`.

The threshold is HARD, deliberately and uniformly. It takes no smoothing width,
because a partial width would make the function non-monotone in the daily range:
smoothing the zero-range case while leaving the general case sharp returns more
exposure at range 0 than at range 1. A hard threshold is also the more faithful
analogue - published sterility criteria use thresholds and ramps, not logistics
- and it is what makes this cell the comparison FLOOR rather than a variant of
the sub-daily paths.

The duration is still smooth in the mean, the range and the threshold through
the `asin`, so the reverse pass is finite; only the saturation points at
`z = +-1` and the degenerate flat day are kinks, which is the same structure as
any clamp in the model.

Because the sub-daily paths smooth their threshold with a 1 C logistic, this
differs from `diurnal_heat_exposure` by that logistic's effect, and the
difference is not small: a day whose maximum sits at the threshold with a 2 C
range accumulates 5.27 h numerically and 0 h here. `test_heat_exposure.jl`
measures that gap instead of tolerating it.
"""
@inline function daily_statistic_exposure_hours(
    mean_temperature::T, amplitude::T, daylength::T, critical_temperature::T,
) where {T <: AbstractFloat}
    daylength <= zero(T) && return zero(T)
    half_range = amplitude * T(0.5)
    # A flat day is all of the window or none of it, on the same hard test as
    # every other day.
    half_range <= eps(T) && return mean_temperature > critical_temperature ?
        daylength : zero(T)
    # Guarded because LLVM speculates `fdiv` out of a branch into a `select`,
    # which would differentiate the division the guard above exists to skip.
    z = guarded_quotient(critical_temperature - mean_temperature, half_range)
    z >= one(T) && return zero(T)
    z <= -one(T) && return daylength
    offset = T(12) * asin(clamp(z, -one(T), one(T))) / T(pi)
    lower = max(T(8) + offset, T(12) - daylength * T(0.5))
    upper = min(T(20) - offset, T(12) + daylength * T(0.5))
    return max(zero(T), upper - lower)
end

"""
    diurnal_heat_exposure(steps, mean_temperature, amplitude, daylength, shape,
                          critical_temperature, width)

Smooth count of daylight hours spent above `critical_temperature`, accumulated
over the same sub-steps the assimilation integral uses.

This is a diagnostic in step 1: it does not affect yield. Its temperature input
becomes organ temperature in step 2 and it begins driving grain number in
step 3. Letting it act on yield while it still reads air temperature would bake
an air-temperature calibration into the reproductive module.
"""
@inline function diurnal_heat_exposure(
    steps::Integer, mean_temperature::T, amplitude::T, daylength::T,
    shape::Integer, critical_temperature::T, width::T,
) where {T <: AbstractFloat}
    daylength <= zero(T) && return zero(T)
    interval = daylength / T(steps)
    total = zero(T)
    for index in 1:steps
        temperature = diurnal_temperature(
            index, steps, mean_temperature, amplitude, daylength, shape,
        )
        total += interval * smooth_exceedance(temperature - critical_temperature, width)
    end
    return total
end
