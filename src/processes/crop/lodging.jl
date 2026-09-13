# Wind lodging: the crop grew, a storm flattened it, and what is on the ground
# is not brought in.
#
# WHY A SECOND WET-YEAR TERM AND NOT A REFINEMENT OF THE FIRST. `excess_water`
# is one of only two mechanisms in this project that improved anything, and its
# own header names lodging among the losses it stands for - on a RAINFALL
# threshold. Whether that is enough is measurable, and it was measured before
# this file was written: over 886 million cropland cell-days, 1981-2016, the
# per-cell year-to-year correlation between the count of high-wind days and the
# count of heavy-rain days has a median of +0.033 (>8 m/s), +0.044 (>10) and
# +0.032 (>12), with |r| > 0.5 on 0.9-1.3% of cells. Wind years and rain years
# are, to a good approximation, different years. The rainfall trigger cannot be
# standing in for this one.
#
# WHY IT CAN WORK WHEN NINE DAMAGE SCALARS DID NOT. `docs/40`: every mechanism
# that lets the plant protect itself damps the anomaly it was built to produce,
# because acclimation is what a plant does to avoid a loss. Shedding leaf area
# saves 24.3% of wheat's drought transpiration and cuts its drought yield loss
# from 41% to 18%. A flattened crop has no such move. Like
# `excess_water_recovery`, this acts inside `harvest_crop!`, at the moment the
# stand is removed, so it cannot feed back into the water balance, the canopy or
# the nitrogen pools - there is no crop left to feed back into.
#
# WHY THE FORCING SUPPORTS IT, WHICH IS NOT OBVIOUS. Berry et al. (2003),
# "Understanding and reducing lodging in cereals", Advances in Agronomy 84, put
# wheat's failure wind speed at 17-20 m/s at ear height. This forcing is a DAILY
# MEAN at 0.5 degrees, and three earlier mechanisms died because the daily mean
# had averaged their event away - 2.8 exposure hours above 35 C in the flowering
# window (`docs/08`), daily rainfall that cannot saturate a profile (`docs/33`),
# and the diurnal course (`docs/26`). Measured on the five gate cells this looked
# like a fourth: maximum daily-mean wind 7.41 m/s over two years. That sample was
# wrong. Over all cropland and 36 years the maximum is 32.79 m/s, 2.98% of
# cell-days exceed 8 m/s and 0.017% exceed 17 m/s. The gate cells are all
# low-wind locations; `tools/diagnostics/wind_ceiling.py` is the measurement that
# corrected it.
#
# WHAT THIS IS NOT. Berry's mechanistic model resolves stem lodging and root
# lodging separately, from stem diameter, wall width, material strength and root
# plate spread. This model carries none of those. What is kept is the structure
# that survives without them: overturning grows with the square of wind speed,
# with the standing load the wind acts on, and with soil wetness, because
# anchorage fails in wet ground. The rate is therefore a calibration of an
# aggregate, not a material property, and the paper should say so.

"""
    lodging!(CFT, crop, wind, soil)

Accumulate today's lodging pressure into the season's exposure.

Runs once per day. Inert before flowering and below the wind threshold, and
inert everywhere when `lodging_rate` is zero, which is the shipped default.
"""
function lodging!(CFT::CFTParameters, crop, wind, soil)
    launch_1D!(
        lodging_kernel!,
        crop_prognostic(crop).phenology.lodging_exposure,
        crop_prognostic(crop).phenology.is_growing,
        crop_phenology_auxiliary(crop).fphu,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.leaf,
        soil_water_auxiliary(soil).relative_content,
        wind,
        CFT,
    )
    return nothing
end

"""
    lodging_pressure_today(wind, wetness, standing_load, fphu, threshold, flowering_start, reference_load)

Today's contribution to the season's lodging exposure, in (m/s)² day.

Three factors, each Berry et al.'s with the material properties this model does
not carry removed:

  - the SQUARE of the wind speed above a threshold, because drag on the canopy
    goes as the square and below the threshold nothing overturns;
  - the standing LOAD the wind acts on, grain plus leaf carbon relative to a
    reference, because an empty canopy has little to overturn;
  - soil WETNESS, because root anchorage fails in wet ground and dry ground
    holds. This is what makes lodging a storm term rather than a wind term.

Zero before flowering: an unfilled ear is light and the stem is still green.
"""
@inline function lodging_pressure_today(
    wind::T, wetness::T, standing_load::T, fphu::T,
    threshold::T, flowering_start::T, reference_load::T,
) where {T <: AbstractFloat}
    fphu >= flowering_start || return zero(T)
    gust = wind - threshold
    gust > zero(T) || return zero(T)
    load = reference_load > zero(T) ?
           min(one(T), max(zero(T), standing_load) / reference_load) : one(T)
    return gust * gust * load * clamp(wetness, zero(T), one(T))
end

"""
    lodging_recovery(season_exposure, tolerance, rate)

The fraction of the standing crop still recoverable after the season's lodging.

One below the tolerance and falling linearly above it, clamped at zero - the
same shape as `excess_water_recovery`, and for the same measured reason. A term
without a tolerance prices a windy CLIMATE rather than a windy YEAR, which is
how a rate large enough to matter in a wet year took 54% of rice yield in an
average one before `heavy_rain_tolerance` existed.
"""
@inline function lodging_recovery(
    season_exposure::T, tolerance::T, rate::T,
) where {T <: AbstractFloat}
    return clamp(one(T) - rate * max(zero(T), season_exposure - tolerance),
                 zero(T), one(T))
end

@kernel inbounds = true function lodging_kernel!(
    lodging_exposure::AbstractVector{T},
    is_growing::AbstractVector{S},
    fphu::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    leaf_carbon::AbstractVector{T},
    relative_content::AbstractMatrix{T},
    wind::AbstractVector{T},
    CFT::CFTParameters,
) where {T <: AbstractFloat, S}
    cell = @index(Global)
    @unpack lodging_wind_threshold, lodging_reference_load, flowering_start = CFT
    if is_growing[cell] != zero(S)
        # Anchorage is a topsoil property: the root plate sits in the first two
        # layers, which hold 95-99% of crop root mass at every gate cell.
        wetness = (relative_content[1, cell] + relative_content[2, cell]) / T(2)
        lodging_exposure[cell] += lodging_pressure_today(
            wind[cell], wetness, storage_carbon[cell] + leaf_carbon[cell],
            fphu[cell], T(lodging_wind_threshold), T(flowering_start),
            T(lodging_reference_load),
        )
    end
end
