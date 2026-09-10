# Drought-driven loss of grain set during flowering.
#
# Water stress already reaches yield twice: through assimilation, via the lambda
# solve and transpiration, and through the harvest index, via
# `compute_harvest_index`'s logistic in season water sufficiency. Neither is the
# mechanism this file adds, and the difference is one word: IRREVERSIBILITY.
#
# The harvest index is recomputed every day from that day's state, so a drought
# that coincides with anthesis and then breaks is forgotten - the model behaves
# as though the rain that came afterwards restored the florets the drought
# aborted. It does not. Water stress at anthesis aborts florets and reduces
# grain NUMBER permanently, which is the same quantity heat sterility reduces
# and the same argument `reproductive_sink.jl` makes for temperature.
#
# So this subtracts from `grain_set_fraction` - the SAME state as heat, not a
# new one, because both stresses abort florets from one pool. Two states
# multiplied would say a drought that halved grain set and a heatwave that
# halved it leave a quarter, which double-counts florets that were already
# gone. Subtraction from a shared, clamped pool is the honest arithmetic, and
# it makes the two drivers directly comparable: the ablation can attribute
# grain-set loss to heat or to water because they are separately switchable
# terms on one state.
#
# The order in which they are applied does not matter. Both only subtract and
# the state is clamped at zero, so `clamp(clamp(x - a) - b)` equals
# `clamp(clamp(x - b) - a)`; `test_water_sterility.jl` asserts it rather than
# trusting it.
#
# WHY THE DAILY SUFFICIENCY AND NOT THE DEFICIT FIELD. `water.sufficiency` is
# today's water scalar in [0, 1] with 1 meaning unstressed, and it is set to 1
# for an absent stand - so a day with no crop reads as no stress, which is the
# sentinel behaviour this mechanism needs. `stress.water_deficit`, despite the
# name, is SEASON-CUMULATIVE sufficiency on 0-100 and is set to ZERO when no
# stand is present. Driving floret abortion from that field would read every
# bare day as total drought, and would integrate an already-integrated
# quantity a second time.
#
# Unlike every other mechanism this project added, this one has no prerequisite:
# `water.sufficiency` is written by `transpiration!` on every day of every
# configuration, so there is no exposure source to enable first.

"""
    water_sterility!(cft, crop)

Reduce `grain_set_fraction` by today's water stress during flowering.

Runs once per day, beside `reproductive_sink!`, on the same contract: the
drivers for today are written and allocation is about to consume the harvest
index that grain set caps.

Inert when `water_sterility_rate` is zero, which is the shipped default until
the rate is bounded against the observational reference the way the heat rates
were - and the bound has to be JOINT with them, because all of these factors
meet in one harvest index.
"""
function water_sterility!(CFT::CFTParameters, crop)
    launch_1D!(
        water_sterility_kernel!,
        crop_prognostic(crop).phenology.grain_set_fraction,
        crop_prognostic(crop).water.sufficiency,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).phenology.is_growing,
        CFT,
    )
    return nothing
end

@kernel inbounds = true function water_sterility_kernel!(
    grain_set_fraction::AbstractVector{T},
    water_sufficiency::AbstractVector{T},
    fphu::AbstractVector{T},
    is_growing::AbstractVector{S},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack flowering_start, flowering_end = CFT
    @unpack water_sterility_sufficiency, water_sterility_rate = CFT

    growing = is_growing[cell] != zero(S)
    # The same raised cosine and the same window as heat sterility: the two
    # stresses abort florets over the same developmental period, so giving them
    # different windows would claim a distinction nothing supports.
    weight = growing ?
        flowering_weight(fphu[cell], T(flowering_start), T(flowering_end)) : zero(T)
    shortfall = max(zero(T), T(water_sterility_sufficiency) - water_sufficiency[cell])
    loss = grain_set_loss(shortfall, weight, T(water_sterility_rate))
    grain_set_fraction[cell] = clamp(grain_set_fraction[cell] - loss, zero(T), one(T))
end
