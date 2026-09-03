"""
fertilizer!(crop, ml, soil, day)

Apply manure/fertilizer inputs and split timing to mineral nitrogen pools.
"""
function fertilizer!(crop,
                     ml::ManagedLand,
                     soil,
                     day;
                     fertilizer::Bool = true,
                     manure::Bool = false,
                     apply_sowing_dose::Bool = true,
                     apply_second_dose::Bool = true,
                     surface_second_manure::Bool = false,
                     reset_inputs::Bool = true,
                     lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (;
        lpjmlparams, fertilizer, manure, apply_sowing_dose,
        apply_second_dose, surface_second_manure, reset_inputs,
    )

    launch_1D!(
        fertilizer_kernel!,
        crop_prognostic(crop).nitrogen.pending_fertilizer,
        crop_calendar_input(crop).sowing_date,
        ml.manure,
        ml.fertilizer,
        crop_prognostic(crop).nitrogen.pending_manure,
        crop_fluxes(crop).nitrogen.prescribed_manure_input,
        crop_fluxes(crop).nitrogen.prescribed_fertilizer_input,
        crop_prognostic(crop).phenology.husum,
        crop_phenology_input(crop).phu,
        soil_nitrogen_prognostic(soil).nitrate,
        soil_nitrogen_prognostic(soil).ammonium,
        soil_carbon_prognostic(soil).litter,
        soil_nitrogen_prognostic(soil).litter,
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
                                    crop_husum::AbstractArray{T},
                                    crop_phu::AbstractArray{T},
                                    soil_NO3::AbstractArray{M},
                                    soil_NH4::AbstractArray{M},
                                    soil_litter_carbon::AbstractArray{M},
                                    soil_litter_nitrogen::AbstractArray{M},
                                    day::Integer,
                                    kernel_params
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)

    @unpack lpjmlparams, fertilizer, manure, apply_sowing_dose,
        apply_second_dose, surface_second_manure, reset_inputs = kernel_params
    @unpack manure_cn, nmanure_nh4_frac, nfert_split_frac, nfert_no3_frac = lpjmlparams

    if reset_inputs
        crop_manure_input[cell] = zero(T)
        crop_fertilizer_input[cell] = zero(T)
    end

    if fertilizer || manure
        fphu = crop_phu[cell] > zero(T) ?
            min(one(T), crop_husum[cell] / crop_phu[cell]) : zero(T)
        if apply_sowing_dose && crop_cal_sdate[cell] == day
            if manure
                manure_input = ml_manure[cell] * nfert_split_frac
                soil_NH4[1, cell] += manure_input * nmanure_nh4_frac
                soil_litter_carbon[2, cell] += manure_input * manure_cn
                soil_litter_nitrogen[2, cell] += manure_input * (1 - nmanure_nh4_frac)
                crop_nmanure[cell] = ml_manure[cell] * (1 - nfert_split_frac)
                crop_manure_input[cell] += manure_input
            end
            if fertilizer
                fertilizer_input = ml_fertilizer[cell] * nfert_split_frac
                soil_NO3[1, cell] += fertilizer_input * nfert_no3_frac
                soil_NH4[1, cell] += fertilizer_input * (1 - nfert_no3_frac)
                crop_nfertilizer[cell] = ml_fertilizer[cell] * (1 - nfert_split_frac)
                crop_fertilizer_input[cell] += fertilizer_input
            end
        end

        if apply_second_dose && fertilizer && fphu > T(0.25) &&
           crop_nfertilizer[cell] > zero(T)
            crop_fertilizer_input[cell] += crop_nfertilizer[cell]
            soil_NO3[1, cell] += crop_nfertilizer[cell] * nfert_no3_frac
            soil_NH4[1, cell] += crop_nfertilizer[cell] * (1 - nfert_no3_frac)
            crop_nfertilizer[cell] = zero(T)
        end

        if apply_second_dose && manure && fphu > T(0.25) &&
           crop_nmanure[cell] > zero(T)
            manure_input = crop_nmanure[cell]
            crop_manure_input[cell] += manure_input
            soil_NH4[1, cell] += manure_input * nmanure_nh4_frac
            litter_pool = surface_second_manure ? 1 : 2
            soil_litter_carbon[litter_pool, cell] += manure_input * manure_cn
            soil_litter_nitrogen[litter_pool, cell] +=
                manure_input * (one(T) - nmanure_nh4_frac)
            crop_nmanure[cell] = zero(T)
        end
    end

end
