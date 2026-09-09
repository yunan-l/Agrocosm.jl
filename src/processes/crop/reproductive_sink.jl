# Heat-driven loss of grain set during flowering.
#
# Every temperature effect in the model reaches the harvest through
# photosynthesis, which means the model cannot produce the signature of
# flowering-period heat damage: yield collapses while biomass barely moves. A
# crop that set no grain still grew. This file adds the missing path.
#
# The design turns on one point. Storage carbon is *recomputed* every day from
# the harvest index, not accumulated, so a sterility factor applied only while
# the crop is hot would be undone the moment the heat passed. Grain that failed
# to set is not recovered by a cool week afterwards, so `grain_set_fraction` is
# a prognostic state that only ever decreases, reset to one at sowing.
#
# Two further choices worth stating, because both are the difference between
# this being a mechanism and being a fudge factor:
#
#   - Exposure is counted as DURATION above a threshold, not as a daily maximum.
#     An hour at 38 C and a day at 38 C are different events and a
#     daily-maximum criterion cannot separate them. This is the same
#     temporal-aggregation argument as the sub-daily assimilation integral,
#     applied to the sink.
#   - The temperature is organ temperature, not air temperature. Feeding air
#     temperature into a sterility threshold would bake an air-temperature
#     calibration into the reproductive module, which is exactly the error this
#     work exists to document.
#
# See docs/05_reproductive_sink_design.md.

# Width of the logistic that replaces the hard "above threshold" indicator, in
# degrees Celsius. `smooth_exceedance` degenerates to the hard test at zero, so
# this is the one knob trading threshold sharpness against a finite reverse
# pass, and both sides of that trade are measured rather than asserted.
#
# It was 1.0, documented as "small enough to behave like a threshold". That
# claim was false and the error was not cosmetic. At width 1.0 the per-hour
# weight is 0.27 at 1 C below the threshold and 0.12 at 2 C below, so a day
# peaking at 34.5 C against a 35 C threshold accumulates 2.95 exposure hours
# having never crossed. Split by whether the day's peak leaf temperature
# actually crossed, that tail supplied 76% of accumulated flowering-window
# exposure at the rice cell, 63% at maize, 44% at soybean and 100% at the
# hot-wheat cell - because those cells have few genuine crossings (rice: 5 days
# against 55) and the tail wins on count. On a uniform sample of days the same
# width contributes only 4%, so the quantity that matters is the width times the
# day distribution, and most cells globally are of the rare-crossing kind.
#
# That biases calibration in the direction that inverts this project's claim: a
# `sterility_rate` fitted against exposure dominated by ordinary warm days, whose
# yields are largely fine, comes out too small, and the mechanism's response to
# real heat events with it. Growing the canopy correctly makes it worse rather
# than better - rice is 22% sub-threshold on the thin, overheated sub-daily
# trajectory and 76% on the calibrated daily one.
#
# 0.3 is where the trade-off turns. The reverse-mode derivative of exposure with
# respect to the daily mean on a hot day, against the weight two degrees below
# the threshold:
#
#   width   weight at -2 C   d(exposure)/d(mean)
#   1.0     0.119            0.842
#   0.5     0.018            0.798
#   0.3     0.0013           0.781     <- here
#   0.2     0.00005          0.742
#   0.1     0                0.476
#   0.05    0                0.092
#
# So 0.3 cuts the sub-threshold weight by about two orders of magnitude while
# keeping 93% of the gradient, and below 0.2 the gradient collapses - which is
# the real constraint that motivated a wide logistic in the first place.
# `docs/07_ablation_framework.md` carries the measurements.
const STERILITY_SMOOTHING_WIDTH = 0.3

"""
    flowering_weight(fphu, start, stop)

Sensitivity of grain set to heat at development stage `fphu`, peaking mid-window
and falling smoothly to zero at either end.

A hard window would put two kinks in the reverse pass and claim a precision
about the timing of anthesis that the heat-unit fraction does not have. Real
susceptibility rises and falls around flowering, so the weight is a raised
cosine over `[start, stop]` and exactly zero outside it.
"""
@inline function flowering_weight(fphu::T, start::T, stop::T) where {T <: AbstractFloat}
    (fphu <= start || fphu >= stop || stop <= start) && return zero(T)
    position = (fphu - start) / (stop - start)
    return T(0.5) * (one(T) - cos(T(2) * T(pi) * position))
end

"""
    grain_set_loss(exposure_hours, weight, rate)

Grain set lost today: linear in accumulated exposure, scaled by developmental
sensitivity.

Linear rather than logistic on purpose. A saturating form needs a shape
parameter that nothing in hand constrains, and it would hide the fact that the
rate itself is uncalibrated. Linear with a floor at zero total set is the
honest starting point and is trivial to replace once validation cells exist.
"""
@inline function grain_set_loss(
    exposure_hours::T, weight::T, rate::T,
) where {T <: AbstractFloat}
    return max(zero(T), rate * exposure_hours * weight)
end

"""
    reproductive_sink!(cft, crop)

Reduce `grain_set_fraction` by today's heat exposure during flowering.

Runs once per day, after the assimilation calls have written
`heat_exposure_hours` at organ temperature. Separate from `photosynthesis!`
because that runs several times a day and accumulating there would count the
same day repeatedly.
"""
function reproductive_sink!(CFT::CFTParameters, crop)
    launch_1D!(
        reproductive_sink_kernel!,
        crop_prognostic(crop).phenology.grain_set_fraction,
        crop_stress_auxiliary(crop).heat_exposure_hours,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.is_growing,
        CFT,
    )
    return nothing
end

@kernel inbounds = true function reproductive_sink_kernel!(
    grain_set_fraction::AbstractVector{T},
    heat_exposure_hours::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack flowering_start, flowering_end, sterility_rate = CFT

    growing = is_growing[cell] != zero(S)
    weight = growing ?
        flowering_weight(fphu[cell], T(flowering_start), T(flowering_end)) : zero(T)
    loss = grain_set_loss(heat_exposure_hours[cell], weight, T(sterility_rate))
    # Monotone by construction: the update only ever subtracts, and the clamp
    # keeps it in [0, 1].
    grain_set_fraction[cell] = clamp(grain_set_fraction[cell] - loss, zero(T), one(T))
end
