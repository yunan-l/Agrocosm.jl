# ifelse is more friendly to GPU parallel computing than idx and @kernel
"""
cultivate!(crop, crop_cal, day)

Handle sowing-day state transitions and activate crop growth state.
"""
function cultivate!(crop::Crop,
                    crop_cal::CropCalendar,
                    ml::ManagedLand,
                    soil::Soil,
                    day::Int;
                    lpjmlparams::LPJmLParams = lpjmlparams,
                    manure = false,
                    apply_prescribed_fertilizer::Bool = true
)

    # if day > 1 && day % 365 == 1
    #     crop_cal.sowing_date = crop_sdate[div(day, 365) + 1, :]
    # end
    # day_ = day % 365 != 0 ? day % 365 : 365
    crop.phenology.harvesting .= ifelse.(crop_cal.sowing_date .== day, false, crop.phenology.harvesting)
    # Update scallback and g_period
    crop_cal.sowing_callback .= ifelse.(crop_cal.sowing_date .== day, 1, crop_cal.sowing_callback)
    crop.phenology.is_growing .= ifelse.(crop_cal.sowing_date .== day, 1, crop.phenology.is_growing)
    crop_cal.sowing_callback .= ifelse.(crop_cal.sowing_date .!= day, 0, crop_cal.sowing_callback)
    fertilizer!(
        crop_cal,
        ml,
        crop,
        soil,
        day;
        enabled = apply_prescribed_fertilizer,
        manure = manure,
        lpjmlparams = lpjmlparams,
    )

    crop.canopy.lai = crop.canopy.lai .* (1 .- crop_cal.sowing_callback) .+ 0.000415f0 .* crop_cal.sowing_callback
    crop.carbon.biomass = crop.carbon.biomass .* (1 .- crop_cal.sowing_callback) .+ 20.0f0 .* crop_cal.sowing_callback
    # init_vegc = device([8.0f0, 0.0113804f0, 0.0f0, 11.9886196f0])
    # crop.carbon.organs = crop.carbon.organs .* (1 .- reshape(crop_cal.sowing_callback, (1, :))) .+ crop.carbon.initial_organs .* reshape(crop_cal.sowing_callback, (1, :))

    # initilization of crop carbon pools and nitrogen pools
    crop.carbon.root = crop.carbon.root .* (1 .- crop_cal.sowing_callback) .+ 8.0f0 .* crop_cal.sowing_callback
    crop.carbon.leaf = crop.carbon.leaf .* (1 .- crop_cal.sowing_callback) .+ 0.0113804f0 .* crop_cal.sowing_callback
    crop.carbon.storage = crop.carbon.storage .* (1 .- crop_cal.sowing_callback) .+ 0.0f0 .* crop_cal.sowing_callback
    crop.carbon.pool = crop.carbon.pool .* (1 .- crop_cal.sowing_callback) .+ 11.9886196f0 .* crop_cal.sowing_callback

    init_nitrogen = 0.7f0 # C:N ratio of seed = 29
    crop.nitrogen.seed_input .= init_nitrogen .* crop_cal.sowing_callback
    crop.nitrogen.total = crop.nitrogen.total .* (1 .- crop_cal.sowing_callback) .+ init_nitrogen * crop_cal.sowing_callback

end
