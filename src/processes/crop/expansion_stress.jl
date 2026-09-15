# WHY A SECOND WATER STRESS. This lineage carries exactly one: `wscal`, which is
# `min(1, supply / demand)` with `supply = emax * f(wr) * (1 - exp(-0.04 rootC))`.
# Two things make it blind over the whole range where a field crop is visibly
# losing yield. `emax` is 8 mm/day for wheat against a demand near 3, and
# `f(wr) = min(1, wr / (1 - p))` is FAO-56's depletion plateau, which is flat
# above `1 - p = 0.45`. At Maricopa the dry arm's root-zone available fraction
# runs at 0.51 through its critical period against the wet arm's 0.93, the field
# lost a fifth of its tillers, and `wscal` reports 1.0000 for BOTH arms.
#
# The plateau is not wrong. It is FAO-56's own statement about TRANSPIRATION: a
# crop transpires at the potential rate until depletion passes `p`. What is wrong
# is using it for everything. Expansive growth is inhibited at water potentials
# well above those that close stomata (Hsiao 1973; Boyer 1970), which is why
# every model that resolves it carries two stresses and not one - CERES/DSSAT's
# `SWDF1` and `SWDF2`, APSIM's `swdef_photo` and `swdef_expansion`, STICS's
# `swfac` and `turfac`. `SWDF2` is `SWDF1` evaluated at 1.5x the demand, so
# expansion is limited at roughly 1.5x the water content transpiration is.
#
# WHAT IT ACTS ON, and why this is grain SET rather than assimilate. Maricopa's
# loss is entirely grain NUMBER: 0.808 and 0.762 of the wet arm across two
# seasons, against a single-grain weight of 1.04 and 1.06 - the dry arm's grains
# were HEAVIER. And the number was lost as tillers, 330 ears against 415 with
# 43.7 grains per ear against 42.9. Tiller number is identical between the arms
# until day 106 and diverges from there, so this is tiller SURVIVAL and not
# tiller production. A reduction in `grain_set_fraction` is exactly that shape:
# fewer grains, each still filling to capacity.

"""
    expansion_stress!(cft, state)

Reduce `grain_set_fraction` by today's shortfall of root-zone available water
below the threshold at which expansive growth is limited, inside the same
flowering window the other grain-set mechanisms use.

`expansion_grain_rate` at zero returns without touching the state, so the shipped
model is unchanged bit for bit.
"""
function expansion_stress!(CFT::CFTParameters, state)
    launch_1D!(
        expansion_stress_kernel!,
        crop_prognostic(state).phenology.grain_set_fraction,
        crop_phenology_auxiliary(state).fphu,
        crop_prognostic(state).phenology.is_growing,
        crop_prognostic(state).water.root_zone_potential,
        CFT,
    )
    return nothing
end

"""
    expansion_stress_loss(weight, threshold, inside_window, rate)

Grain set lost today: linear in the shortfall of the expansive-growth weight
below one, zero outside the window and zero when the crop is unstressed.

THE RATIO AND NOT THE SOIL WATER, which a first attempt got wrong and Braunschweig
caught. Calibrated on root-zone water alone the threshold lands at 0.87, and
Braunschweig's rainfed wheat sits at 0.56 - DRIER than Maricopa's deficit arm at
0.66 - because a humid German season asks for 2.6 mm/day where a Sonoran one asks
for 6.0. On water content alone the mechanism fired on the wrong site and took
that deposit's yield correlation from 0.92 to 0.24. On the ratio the three
separate correctly: 0.89, 1.24 and 1.69.

Linear in the SHORTFALL rather than in a count of stressed days, for the same two
reasons as `anthesis_heat_loss`: a count is not differentiable, and a day at half
the threshold is not a day just below it.
"""
@inline function expansion_stress_loss(
    weight::T, threshold::T, inside_window::Bool, rate::T,
) where {T <: AbstractFloat}
    inside_window || return zero(T)
    isfinite(weight) || return zero(T)
    return max(zero(T), rate * (threshold - weight))
end

@kernel inbounds = true function expansion_stress_kernel!(
    grain_set_fraction::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    root_zone_potential::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack flowering_start, flowering_end = CFT
    @unpack expansion_water_threshold, expansion_grain_rate = CFT

    growing = is_growing[cell] != zero(S)
    # The same rectangular window, open at both ends, as `anthesis_heat_kernel!`,
    # so the two mechanisms damage the same days and compose in any order.
    inside = growing &&
             fphu[cell] > T(flowering_start) &&
             fphu[cell] < T(flowering_end) &&
             T(flowering_end) > T(flowering_start)
    # One where the crop is wetter than the threshold, falling log-linearly to
    # zero at the permanent wilting point.
    weight = expansive_growth_weight(
        root_zone_potential[cell], T(expansion_water_threshold),
    )
    loss = expansion_stress_loss(weight, one(T), inside, T(expansion_grain_rate))
    # Monotone and clamped, like the heat and cold channels, so all three compose
    # in any order.
    grain_set_fraction[cell] = clamp(grain_set_fraction[cell] - loss, zero(T), one(T))
end
