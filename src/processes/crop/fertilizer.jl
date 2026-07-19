"""
fertilizer!(crop_cal, ml, crop, soil, day)

Apply manure/fertilizer inputs and split timing to mineral nitrogen pools.
"""
function fertilizer!(crop_cal::Calendar,
                     ml::Managed_land,
                     crop::Crop,
                     soil::Soil,
                     day;
                     enabled::Bool = true,
                     lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (lpjmlparams = lpjmlparams, enabled = enabled)

    launch_1D!(
        fertilizer_kernel!,
        crop.nfertilizer,
        crop_cal.sdate,
        ml.manure,
        ml.fertilizer,
        crop.nmanure,
        crop.fphu,
        soil.NO3,
        soil.NH4,
        day,
        kernel_params
    )

end


@kernel inbounds = true function fertilizer_kernel!(
                                    crop_nfertilizer::AbstractArray{T},
                                    crop_cal_sdate::AbstractArray{S},
                                    ml_manure::AbstractArray{T},
                                    ml_fertilizer::AbstractArray{T},
                                    crop_nmanure::AbstractArray{T},
                                    crop_fphu::AbstractArray{T},
                                    soil_NO3::AbstractArray{M},
                                    soil_NH4::AbstractArray{M},
                                    day::Integer,
                                    kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack lpjmlparams, enabled = kernel_params
    @unpack nmanure_nh4_frac, nfert_split_frac, nfert_no3_frac = lpjmlparams

    if enabled
        if crop_cal_sdate[cell] == day
            soil_NH4[1, cell] += ml_manure[cell] * nmanure_nh4_frac * nfert_split_frac
            crop_nmanure[cell] = ml_manure[cell] * (1 - nfert_split_frac)

            soil_NO3[1, cell] += ml_fertilizer[cell] * nfert_no3_frac * nfert_split_frac
            soil_NH4[1, cell] += ml_fertilizer[cell] * (1 - nfert_no3_frac) * nfert_split_frac
            crop_nfertilizer[cell] = ml_fertilizer[cell] * (1 - nfert_split_frac)
        end

        if crop_fphu[cell] > T(0.25) && crop_nfertilizer[cell] > zero(T)
            soil_NO3[1, cell] += crop_nfertilizer[cell] * nfert_no3_frac
            soil_NH4[1, cell] += crop_nfertilizer[cell] * (1 - nfert_no3_frac)
            crop_nfertilizer[cell] = zero(T)
        end

        if crop_fphu[cell] > T(0.25) && crop_nmanure[cell] > zero(T)
            soil_NH4[1, cell] += crop_nmanure[cell] * nmanure_nh4_frac
            crop_nmanure[cell] = zero(T)
        end
    end

end
