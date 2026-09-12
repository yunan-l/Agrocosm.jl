"""
temp_stress(CFT, pet, photos, temp; diurnal_range = nothing, steps = 0, shape = 0)

Compute temperature stress scalar used by photosynthesis routines.

With `diurnal_range` supplied the scalar becomes the RADIATION-WEIGHTED mean of
the stress over the day's temperature course instead of the stress at the daily
mean. That is the quantity assimilation actually experiences, and the difference
is large: at a daily mean of 34 C with a 12 C range the retained assimilation is
0.7014 at the mean against 0.3079 over the course for wheat, and 0.9005 against
0.6515 for maize. Rice is unchanged at 1.000 because its `photos.high` of 45 C
puts it off the curve entirely.

This is a change to the RESPONSE FUNCTION, not another damage term, and the
distinction is the measured one. Eleven global arms scored against US county
statistics: nine damage mechanisms, not one of which lowered interannual variance
and the best of which raised correlation by 0.010, against one supply-side fix
that lowered variance and raised correlation on both crops. The sub-daily
assimilation loop, which computes this quantity as a side effect of doing much
more, made soybean substantially worse - r from 0.537 to 0.438 - so this does the
one part that is defensible and none of the rest.

`steps = 0` or a missing range reproduces the daily-mean form bitwise.
"""
function temp_stress(CFT::CFTParameters,
                     pet::PetPar,
                     crop,
                     temp::AbstractArray{T};
                     photoparams::PhotoParams = photoparams,
                     diurnal_range = nothing,
                     diurnal_steps::Integer = 0,
                     diurnal_shape::Integer = DIURNAL_SINUSOID,
) where {T <: AbstractFloat}

    if diurnal_range === nothing || diurnal_steps <= 0
        launch_1D!(
            temp_stress_kernel!,
            crop_photosynthesis_auxiliary(crop).temperature_stress,
            pet.daylength,
            temp,
            CFT,
            photoparams
        )
    else
        launch_1D!(
            diurnal_temp_stress_kernel!,
            crop_photosynthesis_auxiliary(crop).temperature_stress,
            pet.daylength,
            temp,
            diurnal_range,
            CFT,
            photoparams,
            Int32(diurnal_steps),
            Int32(diurnal_shape),
        )
    end

end

@kernel inbounds = true function diurnal_temp_stress_kernel!(
                                     photos_tstress::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     diurnal_range::AbstractArray{T},
                                     CFT::CFTParameters,
                                     photoparams::PhotoParams,
                                     steps::Integer,
                                     shape::Integer,
) where {T <: AbstractFloat}

    cell = @index(Global)

    @unpack path, temp_co2, temp_photos = CFT
    @unpack tmc3, tmc4 = photoparams

    daylength = pet_daylength[cell]
    amplitude = max(zero(T), diurnal_range[cell])
    total = zero(T)
    for step in 1:steps
        # The same course `photosynthesis_subdaily` reconstructs, and the same
        # radiation weights, so the two cannot disagree about what the day looked
        # like - only about how much of the rest of the model runs inside it.
        course = diurnal_temperature(step, steps, temp[cell], amplitude, daylength, shape)
        weight = diurnal_radiation_fraction(step, steps, shape, T)
        total += weight * compute_photosynthesis_temperature_stress(
            daylength, course, path, temp_co2, temp_photos, T(tmc3), T(tmc4),
        )
    end
    photos_tstress[cell] = total
end

"""
    compute_photosynthesis_temperature_stress(daylength, temperature, path, ...)

LPJmL's smooth lower and upper temperature response for C3/C4 assimilation.
The function is scalar and allocation-free so it is safe to call from a CPU or
GPU kernel. `path` uses the existing `1 = C3`, `2 = C4` convention.
"""
@inline function compute_photosynthesis_temperature_stress(
    daylength::T,
    temperature::T,
    path,
    temperature_co2,
    temperature_photosynthesis,
    c3_maximum::T,
    c4_maximum::T,
) where {T <: AbstractFloat}
    # No light, or pathway-specific hard thermal limit: no canopy assimilation.
    (daylength < 0.01 ||
     (path == 1 && temperature > c3_maximum) ||
     (path == 2 && temperature > c4_maximum) ||
     temperature >= temperature_co2.high) && return zero(T)

    lower_slope = T(2 * log(1 / 0.99 - 1)) /
        (temperature_co2.low - temperature_photosynthesis.low)
    # LPJmL fscancftpar.c: midpoint of lower CO2 and photosynthesis limits.
    lower_midpoint = (T(temperature_co2.low) + T(temperature_photosynthesis.low)) * T(0.5)
    upper_slope = T(log(0.99 / 0.01)) /
        (temperature_co2.high - temperature_photosynthesis.high)
    lower_response = 1 / (1 + exp(lower_slope * (lower_midpoint - temperature)))
    upper_response = 1 - 0.01 * exp(
        upper_slope * (temperature - temperature_photosynthesis.high),
    )
    return T(lower_response * upper_response)
end


@kernel inbounds = true function temp_stress_kernel!(
                                     photos_tstress::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     CFT::CFTParameters,
                                     photoparams::PhotoParams
) where {T <: AbstractFloat}

    cell = @index(Global)

    @unpack path, temp_co2, temp_photos = CFT
    @unpack tmc3, tmc4 = photoparams

    photos_tstress[cell] = compute_photosynthesis_temperature_stress(
        pet_daylength[cell], temp[cell], path, temp_co2, temp_photos,
        T(tmc3), T(tmc4),
    )
end
