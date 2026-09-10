# Drought-driven loss of grain filling, after anthesis.
#
# The fourth cell of a 2x2 that the other three mechanisms in this directory
# make: {heat, water} x {grain number, grain weight}.
#
#   reproductive_sink.jl   heat  x number   flowering window, leaf temperature
#   terminal_heat.jl       heat  x weight   filling window, leaf temperature
#   water_sterility.jl     water x number   flowering window, daily sufficiency
#   water_filling.jl       water x weight   filling window, daily sufficiency
#
# The symmetry is not decoration. Measured over the 2015-2016 forcing at the
# five gate cells, the two cells the grain-NUMBER channels cannot touch are both
# reached only through grain WEIGHT, and for the same structural reason in each:
#
#   hot wheat (Punjab)   heat arrives after anthesis. Flowering window holds 2.8
#                        exposure hours above 30 C, the filling window 75.8.
#   wheat (Morocco)      drought arrives after anthesis. Pre-anthesis daily
#                        sufficiency is 1.0 throughout and the flowering window
#                        never drops below 0.539, while the filling window
#                        reaches 0.417 with four days under 0.5.
#
# So a model carrying only grain-number damage is blind at both cells, and which
# stress it is blind to depends on the crop, not on the mechanism.
#
# WHY THIS IS NOT THE EXISTING WATER PATH. `compute_harvest_index` already
# carries a logistic in water sufficiency, so drought already lowers the harvest
# index. That path is driven by `stress.water_deficit`, which is the
# SEASON-CUMULATIVE ratio `sum(min(supply, demand)) / sum(demand)`, and it is
# evaluated at harvest. A mid-season drought does not merely get averaged in
# that quantity, it gets ERASED: every subsequent well-watered day adds the same
# amount to both sums and pulls the ratio back towards one. Measured at the
# Morocco cell, daily sufficiency reaches 0.417 while the season ratio ends at
# 100.0, so the harvest-index penalty there is exactly 100% - no penalty at all
# for a season that contained a drought.
#
# That is this project's own thesis appearing inside the model, on the water
# axis: aggregation destroys the extreme signal, whether the aggregation happens
# in the forcing (a daily mean hiding sub-daily heat) or in the model's own
# bookkeeping (a season ratio hiding a drought). The existing path is left
# untouched - it carries LPJmL's calibration - and this mechanism is the
# non-averaging one placed beside it, so the ablation can measure the
# difference rather than assert it.
#
# See docs/09_water_stress_design.md.

"""
    water_filling!(cft, crop)

Reduce `grain_fill_fraction` by today's water stress during grain filling.

Subtracts from the same state as `terminal_heat!`, because both reduce grain
weight, and for the same reason `water_sterility!` shares a state with
`reproductive_sink!`: two multiplied states would double-count filling capacity
already lost. Order against terminal heat therefore cannot matter.

Runs once per day beside the other three, before allocation consumes the
harvest index. Inert at the shipped rate of zero; the bound has to be taken
jointly with the heat pair, which already had to be scaled to 0.75 of their
separate bounds because all of these factors meet in one harvest index.
"""
function water_filling!(CFT::CFTParameters, crop)
    launch_1D!(
        water_filling_kernel!,
        crop_prognostic(crop).phenology.grain_fill_fraction,
        crop_prognostic(crop).water.sufficiency,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.is_growing,
        CFT,
    )
    return nothing
end

@kernel inbounds = true function water_filling_kernel!(
    grain_fill_fraction::AbstractVector{T},
    water_sufficiency::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack filling_start, filling_end = CFT
    @unpack water_filling_sufficiency, water_filling_rate = CFT

    growing = is_growing[cell] != zero(S)
    weight = growing ?
        flowering_weight(fphu[cell], T(filling_start), T(filling_end)) : zero(T)
    shortfall = max(zero(T), T(water_filling_sufficiency) - water_sufficiency[cell])
    loss = grain_set_loss(shortfall, weight, T(water_filling_rate))
    grain_fill_fraction[cell] = clamp(grain_fill_fraction[cell] - loss, zero(T), one(T))
end
