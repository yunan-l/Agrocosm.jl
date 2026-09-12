# Excess-water damage on an absolute daily RAINFALL threshold, acting on the
# fraction of the standing crop that is recovered at harvest.
#
# WHY A RECOVERY FRACTION AND NOT A STRESS SCALAR. Measured on eleven global
# arms against US county statistics, every damage mechanism this project has
# added raises interannual variance and does not raise correlation - the best of
# nine is +0.010, and the sub-daily loop and canopy energy balance take soybean
# from r = 0.537 to 0.438. The one change that lowers variance AND raises
# correlation is not a damage term at all but a fix to the water supply. So this
# is deliberately not a tenth stress scalar multiplying assimilation: it is a
# pathway the model does not have, and its absence is a SIGN error rather than a
# weak response - in wet disaster cell-years, 12-14% of every crop's disasters,
# the model predicts a good year (+2.0 to +9.8%) where the observation falls
# 22-29%.
#
# WHY IT READS RAINFALL AND NOT SOIL WATER. The obvious mechanism is aeration
# stress on a saturated root zone, which is what APSIM, WOFOST, AquaCrop and
# STICS carry and what Garcia-Vila et al. (2025, Nature Food) review across 21
# field-scale wheat models. It cannot work here, and that is measured rather than
# assumed: under 120 mm/day for five days a 75%-clay column reaches 23% of its
# gravitational pore space and is dry eight days later, because rejected rainfall
# is removed as surface runoff the same day rather than ponding. Anoxia needs
# that pore space 80-90% full.
#
# The pathways that remain are the ones the same review lists beside
# waterlogging - lodging, submergence, sprouting, grain disease, and the field
# access and harvest losses behind the 7.6 Mha left unplanted in the 2019 US
# Midwest. All of them act on the canopy and the grain, and none of them needs a
# saturated soil. A recovery fraction is what they have in common: the crop grew,
# and some of it was not brought in.
#
# WHAT SIZE IT SHOULD BE, OUT OF SAMPLE. Li et al. (2019, Glob. Change Biol.)
# link county USDA crop-insurance indemnities to weather over 1981-2016 - the
# same period, independent of both GDHY and the NASS yields this project scores
# on - and find excessive rainfall costing US maize up to -34% and -17 +/- 3% on
# average, against -37% and -32 +/- 2% for extreme drought. Excessive rainfall is
# the SECOND largest cause of US maize loss after drought, 10 against 18 billion
# dollars over 1989-2016. Two things follow. The mechanism should be able to
# reach a third of yield, not a few percent; and it is not a correction to a
# drought response but a term of comparable size. They also report the loss
# concentrating where drainage is poor, which is the same diagnosis this file's
# second paragraph reaches from the model side.
#
# That band is an out-of-sample CHECK and must stay one: calibrating the rate to
# it would make this a statistical correction wearing a process mechanism's
# clothes.
#
# WHERE THE THRESHOLD COMES FROM. Swept before this file was written, on the
# cell-years that are BOTH a GDHY disaster AND wet. Days above an absolute daily
# rainfall is the best non-circular index for all four crops, and the whole
# season beats every developmental sub-window - which is expected for a term
# that collects lodging early and harvest loss late. `docs/24` has the sweep.

"""
    excess_water!(cft, crop, precipitation)

Reduce `harvest_recovery_fraction` by today's excess of daily rainfall above this
crop's absolute threshold, anywhere in the growing season.
"""
function excess_water!(CFT::CFTParameters, crop, precipitation)
    launch_1D!(
        excess_water_kernel!,
        crop_prognostic(crop).phenology.harvest_recovery_fraction,
        crop_prognostic(crop).phenology.is_growing,
        precipitation,
        CFT,
    )
    return nothing
end

"""
    excess_water_loss(rainfall, threshold, growing, rate)

Recovery lost today: linear in the millimetres by which the day's rainfall
exceeds the threshold, zero below it and zero outside the season.

Linear in the EXCESS rather than counting heavy days, for the reasons
`anthesis_heat_loss` gives: a count is not differentiable, and a 60 mm day is not
a 20.1 mm day.
"""
@inline function excess_water_loss(
    rainfall::T, threshold::T, growing::Bool, rate::T,
) where {T <: AbstractFloat}
    growing || return zero(T)
    return max(zero(T), rate * (rainfall - threshold))
end

@kernel inbounds = true function excess_water_kernel!(
    harvest_recovery_fraction::AbstractVector{T},
    is_growing::AbstractVector{S},
    precipitation::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack heavy_rain_threshold, heavy_rain_rate = CFT
    growing = is_growing[cell] != zero(S)
    loss = excess_water_loss(
        precipitation[cell], T(heavy_rain_threshold), growing, T(heavy_rain_rate),
    )
    # Monotone and clamped, so it composes in any order with anything else that
    # reduces the same state - and irreversible, which is the point: a crop
    # flattened in July is not standing again in September.
    harvest_recovery_fraction[cell] =
        clamp(harvest_recovery_fraction[cell] - loss, zero(T), one(T))
end
