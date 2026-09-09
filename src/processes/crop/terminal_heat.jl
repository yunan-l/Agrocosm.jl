# Heat damage to grain FILLING, after anthesis.
#
# The reproductive sink in `reproductive_sink.jl` fixes how many grains are set,
# around flowering. This file adds the other half of heat damage to yield: how
# much each grain fills, afterwards. They are not variants of one mechanism and
# a model carrying only the first cannot express terminal heat at all, which is
# the canonical wheat failure - a crop that set a full complement of grains and
# then filled them poorly because the last month was hot.
#
# Two measurements at the hot-wheat gate cell (Indian Punjab, the textbook
# terminal-heat case) forced this to be separate rather than a retuning of the
# sink, and both are in `docs/08_terminal_heat_design.md`:
#
#   - WRONG WINDOW. The heat there arrives after anthesis. The flowering window
#     accumulates 2.8 leaf-temperature exposure hours above 30 C; the filling
#     window accumulates 75.8.
#   - WRONG THRESHOLD. Peak leaf temperature in that filling window is 34.5 C,
#     so exposure above the 35 C sterility threshold is 0.4 h. The sink is inert
#     there at every rate up to 0.1 - measured, not assumed - and no rate can
#     fix a threshold that is never crossed. Grain filling is impaired well
#     below sterility, so `filling_temperature` is lower and
#     `filling_exposure_hours` is accumulated against it separately.
#
# Design follows the sink's, deliberately, so the two can be compared:
#
#   - `grain_fill_fraction` is PROGNOSTIC and only ever decreases, reset to one
#     at sowing. Storage carbon is recomputed daily from the harvest index, so a
#     factor applied only while the crop is hot would be undone by the next cool
#     week - and starch not deposited is not deposited later.
#   - Exposure is DURATION above a threshold, not a daily maximum.
#   - The temperature is organ temperature, for the same reason: an
#     air-temperature calibration baked into a reproductive module is the error
#     this project exists to document.
#
# What this does NOT represent, to be stated in the paper rather than implied
# away: filling damage is lumped into one multiplier on the harvest index, so
# the model does not distinguish a shortened filling DURATION from a suppressed
# filling RATE, and it has no stem-reserve remobilisation to fail. The
# developmental weight is the same raised cosine the flowering window uses,
# which puts peak sensitivity mid-filling; the physiological evidence is closer
# to flat-then-declining, so the shape is a stand-in chosen for continuity with
# the sink and for a kink-free reverse pass, not from filling data.

"""
    terminal_heat!(cft, crop)

Reduce `grain_fill_fraction` by today's heat exposure during grain filling.

Runs once per day, immediately after `reproductive_sink!`, on the same contract:
the exposure fields have been written for today, and the harvest index that
allocation is about to use has to already reflect both losses.

Inert when `filling_rate` is zero, which is the shipped default - the window
carries far more exposure than the flowering window at every gate cell, so a
nonzero rate is a large intervention and is bounded per crop before it is set.
"""
function terminal_heat!(CFT::CFTParameters, crop)
    launch_1D!(
        terminal_heat_kernel!,
        crop_prognostic(crop).phenology.grain_fill_fraction,
        crop_stress_auxiliary(crop).filling_exposure_hours,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.is_growing,
        CFT,
    )
    return nothing
end

@kernel inbounds = true function terminal_heat_kernel!(
    grain_fill_fraction::AbstractVector{T},
    filling_exposure_hours::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack filling_start, filling_end, filling_rate = CFT

    growing = is_growing[cell] != zero(S)
    # The same raised cosine as the flowering window, over the filling window:
    # exactly zero at both ends, so there is no kink where the windows meet.
    weight = growing ?
        flowering_weight(fphu[cell], T(filling_start), T(filling_end)) : zero(T)
    loss = grain_set_loss(filling_exposure_hours[cell], weight, T(filling_rate))
    # Monotone by construction, like grain set: the update only subtracts.
    grain_fill_fraction[cell] = clamp(grain_fill_fraction[cell] - loss, zero(T), one(T))
end
