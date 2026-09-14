# Loss of stand to heat between sowing and a completed canopy.
#
# WHY THIS EXISTS AT ALL. Every temperature effect this model carries acts at or
# after flowering: `anthesis_heat!`, `cold_sterility!`, `terminal_heat!`,
# `water_sterility!`, and the temperature limb on assimilation. Nothing acts
# between sowing and canopy closure, so a crop in this lineage CANNOT fail to
# establish, and therefore cannot produce a zero-yield season from heat before
# flowering however hot the season is.
#
# The Hot Serial Cereal experiment sowed wheat on twelve dates across the year at
# Maricopa. Three sowings produced NO GRAIN AT ALL - no ears, no anthesis date,
# nothing to measure - and they are exactly the three whose first thirty days
# carried 25 to 30 days above 35 C, against 0 to 12 for every sowing that
# yielded. The model returns 0.77, 0.85 and 4.15 t/ha for them.
#
# CARDINAL TEMPERATURE, NOT A FITTED THRESHOLD. Porter and Gawith (1999,
# Eur. J. Agron.) review wheat's cardinal temperatures and give germination a
# maximum of 35 C. That is the ceiling below.
#
# WHAT IT ACTS ON. Stand, expressed as a multiplier on POTENTIAL leaf area: a
# crop that has lost half its plants builds half the canopy, and one that has
# lost all of them builds none. That is the physical quantity, and it is the
# reason this cannot be repaired later - the plants are dead. `stand_fraction` is
# therefore prognostic and monotone, reset to one at sowing, exactly as
# `grain_set_fraction` is.
#
# SHIPS AT ZERO RATE, which is the inherited model bitwise.

"""
    establishment_loss(temperature, ceiling, rate, width)

The share of the stand a day above `ceiling` destroys.

Smoothed with the same logistic as the reproductive sink so the reverse pass
stays finite, and at `rate = 0` it returns zero whatever the day was.
"""
@inline function establishment_loss(temperature::T, ceiling::T, rate::T,
                                    width::T) where {T <: AbstractFloat}
    rate > zero(T) || return zero(T)
    return rate * smooth_exceedance(temperature - ceiling, width)
end

"""
    establishment!(CFT, crop)

Accumulate heat loss of stand while the canopy is still being built.

The window is a fixed number of DAYS after sowing. Closing it on thermal time
instead closes it fastest in the hottest season - 5 days against 13 across these
sowings - which lets the crop escape the stress the window represents.
"""
function establishment!(CFT::CFTParameters, crop, temperature::AbstractArray,
                        diurnal_range::AbstractArray)
    launch_1D!(
        establishment_kernel!,
        crop_prognostic(crop).phenology.stand_fraction,
        crop_prognostic(crop).phenology.growing_days,
        crop_prognostic(crop).phenology.is_growing,
        temperature,
        diurnal_range,
        (ceiling = CFT.establishment_heat_ceiling,
         rate = CFT.establishment_loss_rate,
         window_days = CFT.establishment_days,
         width = STERILITY_SMOOTHING_WIDTH),
    )
    return nothing
end

@kernel inbounds = true function establishment_kernel!(
    stand_fraction::AbstractArray{T},
    growing_days::AbstractArray{S},
    is_growing::AbstractArray{S},
    temperature::AbstractArray{T},
    diurnal_range::AbstractArray{T},
    parameters,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if is_growing[cell] == one(S) && growing_days[cell] <= S(parameters.window_days)
        # The DAILY MAXIMUM, not the mean. Porter and Gawith's ceiling is a
        # process limit, and a seedbed reaches it in the afternoon: across the
        # twelve Hot Serial Cereal sowings the daily maximum separates the three
        # total failures from every harvested treatment (25-30 days above 35 C in
        # the first thirty against 0-12), while the daily mean of the failures,
        # 28.8 to 32.5 C, never crosses it at all.
        daily_maximum = temperature[cell] + max(diurnal_range[cell], zero(T)) / T(2)
        loss = establishment_loss(daily_maximum, T(parameters.ceiling),
                                  T(parameters.rate), T(parameters.width))
        stand_fraction[cell] = max(stand_fraction[cell] - loss, zero(T))
    end
end
