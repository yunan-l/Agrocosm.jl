# Cold sterility on an absolute daily-minimum AIR temperature, the mirror of
# `anthesis_heat!`.
#
# WHY THIS EXISTS AT ALL. The model has no cold damage of any kind. It carries
# vernalization, which is a REQUIREMENT for development rather than an injury,
# and a low-temperature limb on assimilation - and that limb is evaluated on the
# DAILY MEAN, which is the same defect the heat side has. A rice booting night at
# 12 C, well below every published cold-sterility threshold, arrives at the
# response function as a daily mean near 19 C and is answered with a 1.9%
# assimilation penalty (`temp_stress.jl:52`, rice `temp_co2.low = 6`,
# `temp_photos.low = 20`). Worse, that penalty is on SOURCE, which the season can
# make up, while cold at microsporogenesis destroys SINK, which it cannot.
#
# This is a gap in the lineage, not in this implementation, and it is a gap the
# field has documented: frost enters most crop models through winterkill
# functions acting on the stand, and of the models reviewed by Barlow et al.
# (2015, Field Crops Res.) only STICS reduces GRAIN NUMBER around anthesis.
# ORYZA v3 is the other exception, simulating chilling sterility from cooling
# degree-days. Neither route is in LPJmL and so neither is in Agrocosm.
#
# WHAT THE MEASUREMENT SAYS, AND WHAT IT DOES NOT. Swept exactly as the heat
# thresholds were - days below an absolute daily minimum, inside a developmental
# window, z-scored per cell, averaged over the cell-years GDHY calls disasters.
# The cold signal is REAL and SMALL: maize peaks at +0.363 (booting, <2 C)
# against its heat peak of +1.076, and rice's best cold number, +0.394, rests on
# 805 of 15,704 disaster cell-years. It is a regional mechanism carried by
# high-latitude cells, and it cannot be expected to close the amplitude gap.
# The honest statement is that this closes a STRUCTURAL hole - the model could
# not represent cold injury at all - not that it closes the amplitude gap.
# `docs/19` records the full sweep including the thresholds it failed to
# confirm.

"""
    cold_sterility!(cft, crop, temperature, diurnal_range)

Reduce `grain_set_fraction` by today's shortfall of DAILY MINIMUM air temperature
below this crop's absolute cold-sterility threshold, inside the microsporogenesis
window.

Requires `diurnal_range` in the forcing, for the same reason `anthesis_heat!`
does: a daily minimum cannot be recovered from a mean.
"""
function cold_sterility!(CFT::CFTParameters, crop, temperature, diurnal_range)
    launch_1D!(
        cold_sterility_kernel!,
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
    cold_sterility_loss(daily_minimum, threshold, inside_window, rate)

Grain set lost today: linear in the degrees by which the daily minimum falls
below the threshold, zero outside the window and zero above it.

Accumulated over the window by repeated subtraction, this is a COOLING
DEGREE-DAY response, the form ORYZA v3 and SIMRIW use for chilling sterility -
but on the daily minimum rather than the daily mean, because the minimum is
where the published thresholds are stated and because the mean is the quantity
the model's existing low-temperature limb already sees.
"""
@inline function cold_sterility_loss(
    daily_minimum::T, threshold::T, inside_window::Bool, rate::T,
) where {T <: AbstractFloat}
    inside_window || return zero(T)
    return max(zero(T), rate * (threshold - daily_minimum))
end

@kernel inbounds = true function cold_sterility_kernel!(
    grain_set_fraction::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    temperature::AbstractVector{T},
    diurnal_range::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack cold_start, cold_end, cold_night_temperature, cold_night_rate = CFT

    growing = is_growing[cell] != zero(S)
    # A window of its own, NOT `flowering_start`/`flowering_end`. Cold sterility
    # acts at microsporogenesis, which precedes anthesis by roughly eleven days
    # in rice; heat sterility acts at anthesis itself. Sharing one window would
    # be wrong physiology even where the two happen to overlap. Rectangular and
    # open at both ends, matching `anthesis_heat!`.
    inside = growing &&
             fphu[cell] > T(cold_start) &&
             fphu[cell] < T(cold_end) &&
             T(cold_end) > T(cold_start)
    # Daily minimum from the two statistics the forcing carries. `diurnal_range`
    # is tasmax - tasmin, so half of it below the mean is the minimum.
    daily_minimum = temperature[cell] - max(zero(T), diurnal_range[cell]) * T(0.5)
    loss = cold_sterility_loss(
        daily_minimum, T(cold_night_temperature), inside, T(cold_night_rate),
    )
    # Monotone and clamped, so this composes with the four mechanisms that share
    # `grain_set_fraction` in any order.
    grain_set_fraction[cell] = clamp(grain_set_fraction[cell] - loss, zero(T), one(T))
end
