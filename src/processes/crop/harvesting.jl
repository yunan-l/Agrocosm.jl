"""
harvest_crop!(crop_cal, crop, soil, day)

Handle harvest-day biomass removal, residue transfer, and crop state reset.
"""
function harvest_crop!(crop_cal::CropCalendar,
                       crop::Crop,
                       soil::Soil,
                       output::Output,
                       residue_frac::AbstractArray{T},
                       day::Int
) where {T <: AbstractFloat}

    # update hcallback and g_period
    crop_cal.harvest_date .= ifelse.((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== true), day, crop_cal.harvest_date)
    crop_cal.harvest_callback .= ifelse.((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== true), 1, 0)
    crop.phenology.is_growing .= ifelse.((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== true), 0, crop.phenology.is_growing)
    crop.nitrogen.harvest_export .=
        (crop.nitrogen.storage .+
         (crop.nitrogen.leaf .+ crop.nitrogen.pool) .* (one(T) .- residue_frac)) .*
        crop_cal.harvest_callback
    # Update crop variables
    crop.carbon.yield .= ifelse.(((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== true)), crop.carbon.storage, crop.carbon.yield)

    soil.carbon.input[SURFACE_LITTER, :] .= ((crop.carbon.leaf .+ crop.carbon.pool) .* residue_frac) .* crop_cal.harvest_callback
    soil.carbon.input[INCORPORATED_LITTER, :] .= zero(T)
    soil.carbon.input[ROOT_LITTER, :] .= crop.carbon.root .* crop_cal.harvest_callback

    soil.nitrogen.input[SURFACE_LITTER, :] .= ((crop.nitrogen.leaf .+ crop.nitrogen.pool) .* residue_frac) .* crop_cal.harvest_callback
    soil.nitrogen.input[INCORPORATED_LITTER, :] .= zero(T)
    soil.nitrogen.input[ROOT_LITTER, :] .= crop.nitrogen.root .* crop_cal.harvest_callback

    # Do not clear plant N here. harvest_crop! has already set is_growing = 0,
    # so the later same-day crop_nitrogen! call clears total N in
    # nuptake_crop! and all organ N pools in allocate_crop_nitrogen!. Keeping
    # that operation in the kernels avoids five redundant GPU broadcasts.
    # crop.nitrogen.total .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.leaf .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.root .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.storage .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.pool .*= one(T) .- crop_cal.harvest_callback

    # soil.carbon.input .= vcat(reshape((crop.carbon.leaf .+ crop.carbon.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.carbon.root, (1, :))) .* reshape(crop_cal.harvest_callback, (1, :))
    # soil.nitrogen.input .= vcat(reshape((crop.nitrogen.leaf .+ crop.nitrogen.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.nitrogen.root, (1, :))) .* reshape(crop_cal.harvest_callback, (1, :))
    # idx = ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== true)) .| ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== false)) .| ((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== false))
    # crop_cal.harvest_callback[idx] .= 0
    # crop_cal.harvest_callback .= ifelse.(((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== true)) .| ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== false)) .| ((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== false)), 0, crop_cal.harvest_callback)

    # update harvesting variables
    output.crop.growing_mask = vcat(output.crop.growing_mask, reshape(crop.phenology.is_growing, (1, :)))
    output.calendar.harvesting_mask = vcat(output.calendar.harvesting_mask, reshape(crop_cal.harvest_callback, (1, :)))
    output.crop.storage_carbon = vcat(output.crop.storage_carbon, reshape(crop.carbon.storage, (1, :)))
    output.calendar.sowing_callback = vcat(output.calendar.sowing_callback, reshape(crop_cal.sowing_callback, (1, :)))
    output.calendar.harvest_callback = vcat(output.calendar.harvest_callback, reshape(crop_cal.harvest_callback, (1, :)))
    output.crop.water_deficit = vcat(output.crop.water_deficit, reshape(mean((soil.water.storage ./ soil.properties.layer_depth)[1:3, :], dims = 1), (1, :)))
    if day == 365
        output.calendar.harvest_date = vcat(output.calendar.harvest_date, reshape(crop_cal.harvest_date, (1, :)))
        crop_cal.harvesting_year .= ifelse.(crop.carbon.yield .!= 0.0f0, 1, 0)
        output.crop.yield = vcat(output.crop.yield, reshape(max.(crop.carbon.yield, 0.0f0), (1, :)))
        crop.carbon.yield .= 0.0f0
        crop_cal.harvest_date .= 0
        output.calendar.harvesting_year = vcat(output.calendar.harvesting_year, reshape(crop_cal.harvesting_year, (1, :)))
    end
    # crop.carbon.organs = crop.carbon.organs .* (1 .- reshape(crop_cal.harvest_callback, (1, :)))
    # crop.carbon.root = crop.carbon.root .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.leaf = crop.carbon.leaf .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.storage = crop.carbon.storage .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.pool = crop.carbon.pool .* (1 .- crop_cal.harvest_callback)

end
