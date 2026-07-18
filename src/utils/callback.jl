"""
update_lit_winter_wheat!(soil, litch, litnh, crop_wtype, hdate, crop_cal_hcallback, day)

Apply winter-wheat harvest callback logic for litter pool updates.
"""
function update_lit_winter_wheat!(soil::Soil,
                                  litch::AbstractArray{M},
                                  litnh::AbstractArray{M},
                                  crop_wtype::AbstractArray{B},
                                  hdate::AbstractArray{S},
                                  crop_cal_hcallback::AbstractArray{S},
                                  day::Int
) where {M <: AbstractFloat, B <: Bool, S <: Integer}

    hdate_callback = copy(crop_cal_hcallback)

    Zygote.ignore() do
        launch_1D!(update_lit_winter_wheat_kernel!,
                   hdate_callback,
                   crop_wtype,
                   hdate,
                   day)
    end

    soil.litc = soil.litc .* (1 .- reshape(hdate_callback, (1, :))) + litch .* reshape(hdate_callback, (1, :)) 
    soil.litn = soil.litn .* (1 .- reshape(hdate_callback, (1, :))) + litnh .* reshape(hdate_callback, (1, :)) 

end


@kernel inbounds = true function update_lit_winter_wheat_kernel!(
                                                 hdate_callback::AbstractArray{S},
                                                 crop_wtype::AbstractArray{B},
                                                 hdate::AbstractArray{S},
                                                 day::Int

) where {B <: Bool, S <: Integer}

    cell = @index(Global)

    hdate_callback[cell] = zero(S)

    if day < 365 && (crop_wtype[cell] == true) && (hdate[cell] == day)
        hdate_callback[cell] = one(S)
    end

end
# Only used for winter wheat training!!!