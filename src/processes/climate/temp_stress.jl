"""
temp_stress(PFT, pet, photos, temp)

Compute temperature stress scalar used by photosynthesis routines.
"""
function temp_stress(PFT::PftParameters,
                     pet::PetPar,
                     crop::Crop,
                     temp::AbstractArray{T};
                     photoparams::PhotoParams = photoparams
) where {T <: AbstractFloat}

    launch_1D!(
        temp_stress_kernel!,
        crop.auxiliary.photosynthesis.temperature_stress,
        pet.daylength,
        temp,
        PFT,
        photoparams
    )

end


@kernel inbounds = true function temp_stress_kernel!(
                                     photos_tstress::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     PFT::PftParameters,
                                     photoparams::PhotoParams
) where {T <: AbstractFloat}

    cell = @index(Global)

    @unpack path, temp_co2, temp_photos = PFT
    @unpack tmc3, tmc4 = photoparams

    k1 = T(2 * log(1 / 0.99 - 1)) / (temp_co2.low - temp_photos.low)
    k2 = temp_co2.low + temp_photos.low * T(0.5)
    k3 = T(log(0.99 / 0.01)) / (temp_co2.high - temp_photos.high)

    if pet_daylength[cell] < 0.01 || (path == 1 && temp[cell] > tmc3) || (path == 2 && temp[cell] > tmc4) # path == 1 : C3; path == 2 : C4
        photos_tstress[cell] = zero(T)
    else
        if temp[cell] < temp_co2.high
            low = 1 / (1 + exp(k1 * (k2 - temp[cell])))
            high = 1 - 0.01 .* exp(k3 * (temp[cell] - temp_photos.high))
            photos_tstress[cell] = T(low * high)
        else
            photos_tstress[cell] = zero(T)
        end
    end
end
