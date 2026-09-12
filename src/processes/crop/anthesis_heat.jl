# Anthesis heat sterility on an absolute daily-maximum AIR temperature.
#
# WHY A FOURTH HEAT PATH WHEN THREE ALREADY EXIST. `reproductive_sink!` reads
# accumulated exposure HOURS at ORGAN temperature above a single 35 C threshold
# shared by every crop. Measured globally against GDHY, that combination does not
# work, and each of its three choices is wrong for a separate reason:
#
#   * organ temperature is COOLER than air when the canopy transpires, so adding
#     it makes the model LESS responsive in disaster years - soybean moves from
#     -0.623 to -0.482 and rice from -0.196 to -0.167 on the ablation ladder;
#   * accumulated hours re-average inside the day the threshold was meant to
#     pick out, and on the same ladder `tmax_sink`, which reads a daily
#     statistic directly, beats both `subdaily` and `organ_temperature`;
#   * one threshold for four crops cannot be right when the measured optima are
#     36, 38, 38 and 38 C.
#
# WHERE THE THRESHOLDS COME FROM. They were measured before this file was
# written, not chosen and then bounded: a sweep of exceedance-day counts over the
# cell-years GDHY calls disasters, per crop, per developmental window, on the
# same `tasmax` the model reads. All four crops peak in the FLOWERING window, and
# maize peaks at 38 C against the 37.9 C that controlled-environment pollen work
# reports. See `docs/18_heat_day_assimilation_design.md` for the full sweep and
# for wheat's, which is flat and is the one number here not to trust.
#
# WHY IT SUBTRACTS FROM GRAIN SET AND NOT FROM ASSIMILATION. Sterile pollen
# removes SINK capacity, not source. `docs/18` first proposed an assimilation
# penalty on the grounds that the harvest index is the binding constraint on only
# 38.5% of growing days in wheat and rice, so a multiplier on it is often inert.
# That is an argument from architecture, not from physiology, and it is also
# wrong on its own terms: `compute_storage_carbon` takes `min(sink, source)`, so
# a sink reduction bites whenever it goes BELOW the source - it is small rates
# that fail to bite, not the pathway. The risk this leaves is real and is stated
# rather than designed around: if wheat and rice do not respond at a defensible
# rate, the mass cap is the reason and the assimilation route is the fallback.

"""
    anthesis_heat!(cft, crop, temperature, diurnal_range)

Reduce `grain_set_fraction` by today's excess of DAILY MAXIMUM air temperature
above this crop's absolute sterility threshold, inside the flowering window.

Requires `diurnal_range` in the forcing, like the daily-statistic exposure path:
a plain daily run carries no such field, and the daily maximum cannot be
recovered from a mean alone.
"""
function anthesis_heat!(CFT::CFTParameters, crop, temperature, diurnal_range)
    launch_1D!(
        anthesis_heat_kernel!,
        crop_prognostic(crop).phenology.grain_set_fraction,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.is_growing,
        temperature,
        diurnal_range,
        CFT,
    )
    return nothing
end

"""
    anthesis_heat_loss(daily_maximum, threshold, inside_window, rate)

Grain set lost today: linear in the degrees by which the daily maximum exceeds
the threshold, zero outside the window and zero below it.

Linear in the EXCESS rather than counting exceedance days, because a count is not
differentiable - which this project cannot afford - and because a day at 42 C is
not a day at 38.1 C, which a count cannot express.
"""
@inline function anthesis_heat_loss(
    daily_maximum::T, threshold::T, inside_window::Bool, rate::T,
) where {T <: AbstractFloat}
    inside_window || return zero(T)
    return max(zero(T), rate * (daily_maximum - threshold))
end

@kernel inbounds = true function anthesis_heat_kernel!(
    grain_set_fraction::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    temperature::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack flowering_start, flowering_end = CFT
    @unpack heat_day_temperature, heat_day_rate = CFT

    growing = is_growing[cell] != zero(S)
    # A RECTANGULAR window, not the raised cosine the other reproductive
    # mechanisms use. That cosine was chosen for a kink-free reverse pass on a
    # mechanism whose window nothing had measured; this window WAS measured, and
    # its edges are where the signal stops rather than where a smooth weight
    # happens to vanish. Open at both ends, matching `flowering_weight`.
    inside = growing &&
             fphu[cell] > T(flowering_start) &&
             fphu[cell] < T(flowering_end) &&
             T(flowering_end) > T(flowering_start)
    # Daily maximum from the two statistics the forcing carries. `diurnal_range`
    # is tasmax - tasmin, so half of it above the mean is the maximum.
    daily_maximum = temperature[cell] + max(zero(T), diurnal_range[cell]) * T(0.5)
    loss = anthesis_heat_loss(
        daily_maximum, T(heat_day_temperature), inside, T(heat_day_rate),
    )
    # Monotone by construction and clamped, so this composes with the other three
    # mechanisms in any order: clamp(clamp(x - a) - b) == clamp(clamp(x - b) - a).
    grain_set_fraction[cell] = clamp(grain_set_fraction[cell] - loss, zero(T), one(T))
end
