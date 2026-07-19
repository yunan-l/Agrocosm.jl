"""
fertilizer!(crop_cal, ml, crop, soil, day)

Apply manure/fertilizer inputs and split timing to mineral nitrogen pools.
"""
function fertilizer!(crop_cal::CropCalendar,
                     ml::ManagedLand,
                     crop::Crop,
                     soil::Soil,
                     day;
                     enabled::Bool = true,
                     manure::Bool = false,
                     lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (; lpjmlparams, enabled, manure)

    launch_1D!(
        fertilizer_kernel!,
        crop.nitrogen.pending_fertilizer,
        crop_cal.sowing_date,
        ml.manure,
        ml.fertilizer,
        crop.nitrogen.pending_manure,
        crop.nitrogen.prescribed_manure_input,
        crop.nitrogen.prescribed_fertilizer_input,
        crop.phenology.fphu,
        soil.nitrogen.nitrate,
        soil.nitrogen.ammonium,
        soil.carbon.litter,
        soil.nitrogen.litter,
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
                                    crop_manure_input::AbstractArray{T},
                                    crop_fertilizer_input::AbstractArray{T},
                                    crop_fphu::AbstractArray{T},
                                    soil_NO3::AbstractArray{M},
                                    soil_NH4::AbstractArray{M},
                                    soil_litter_carbon::AbstractArray{M},
                                    soil_litter_nitrogen::AbstractArray{M},
                                    day::Integer,
                                    kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)

    @unpack lpjmlparams, enabled, manure = kernel_params
    @unpack manure_cn, nmanure_nh4_frac, nfert_split_frac, nfert_no3_frac = lpjmlparams

    crop_manure_input[cell] = zero(T)
    crop_fertilizer_input[cell] = zero(T)

    if enabled
        if crop_cal_sdate[cell] == day
            fertilizer_input = ml_fertilizer[cell] * nfert_split_frac
            if manure
                manure_input = ml_manure[cell] * nfert_split_frac
                soil_NH4[1, cell] += manure_input * nmanure_nh4_frac
                soil_litter_carbon[2, cell] += manure_input * manure_cn
                soil_litter_nitrogen[2, cell] += manure_input * (1 - nmanure_nh4_frac)
                crop_nmanure[cell] = ml_manure[cell] * (1 - nfert_split_frac)
                crop_manure_input[cell] += manure_input
            end

            soil_NO3[1, cell] += ml_fertilizer[cell] * nfert_no3_frac * nfert_split_frac
            soil_NH4[1, cell] += ml_fertilizer[cell] * (1 - nfert_no3_frac) * nfert_split_frac
            crop_nfertilizer[cell] = ml_fertilizer[cell] * (1 - nfert_split_frac)
            crop_fertilizer_input[cell] += fertilizer_input
        end

        if crop_fphu[cell] > T(0.25) && crop_nfertilizer[cell] > zero(T)
            crop_fertilizer_input[cell] += crop_nfertilizer[cell]
            soil_NO3[1, cell] += crop_nfertilizer[cell] * nfert_no3_frac
            soil_NH4[1, cell] += crop_nfertilizer[cell] * (1 - nfert_no3_frac)
            crop_nfertilizer[cell] = zero(T)
        end

        if manure && crop_fphu[cell] > T(0.25) && crop_nmanure[cell] > zero(T)
            manure_input = crop_nmanure[cell]
            crop_manure_input[cell] += manure_input
            soil_NH4[1, cell] += manure_input * nmanure_nh4_frac
            soil_litter_carbon[2, cell] += manure_input * manure_cn
            soil_litter_nitrogen[2, cell] += manure_input * (1 - nmanure_nh4_frac)
            crop_nmanure[cell] = zero(T)
        end
    end

end
