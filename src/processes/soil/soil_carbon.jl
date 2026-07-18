"""
soil_carbon!(crop_cal, crop, soil)

Update litter and soil carbon pools and heterotrophic respiration terms.
"""
function soil_carbon!(crop_cal::Calendar,
                      soil::Soil;
                      lpjmlparams::LPJmLParams = lpjmlparams,
                      soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack atmfrac, k_soil10 = lpjmlparams
    @unpack e0, intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # soil decomposition response
    soil_decomp_response!(soil)

    # compute soil carbon: litter carbon and soil carbon
    # soil.decom_litc = (1.0f0 .- exp.(-soil.response_litc / 100)) .* soil.litc

    # Litter decomposition is represented as three aggregated litter pools.
    # We use top-layer response (LPJmL uses top/root layer litter environments).
    soil.decom_litc = (1.0f0 .- exp.(-soil.response_litc .* soil.decom_lit_response)) .* soil.litc
    soil.litc = soil.litc  - soil.decom_litc

    # using 'callback' to adjust litter carbon due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    update_litc_tillage!(soil, crop_cal)
    
    # soil.decom_fastc = (1.0f0 .- exp.(-soil.response_fastc .* response / 50)) .* soil.fastc
    soil.decom_fastc = (1.0f0 .- exp.(-k_soil10.fast .* soil.decom_response)) .* soil.fastc
    soil.fastc = soil.fastc + soil.c_shift_fast .* sum(soil.decom_litc, dims = 1) - soil.decom_fastc
    
    # soil.decom_slowc = (1.0f0 .- exp.(-soil.response_slowc .* response / 10)) .* soil.slowc
    soil.decom_slowc = (1.0f0 .- exp.(-k_soil10.slow .* soil.decom_response)) .* soil.slowc
    soil.slowc = soil.slowc + soil.c_shift_slow .* sum(soil.decom_litc, dims = 1) - soil.decom_slowc

    soil.rh = vec(sum(soil.decom_litc, dims = 1) * atmfrac .+ sum(soil.decom_fastc, dims = 1) .+ sum(soil.decom_slowc, dims = 1))
    
end


"""
update_litc_tillage!(soil, crop_cal)

Apply tillage/harvest crop carbon to litter carbon pools.
"""
function update_litc_tillage!(soil::Soil,
                              crop_cal::Calendar
)

    soil.litc = soil.litc .* (1 .- reshape(crop_cal.scallback, (1, :))) .* (1 .- reshape(crop_cal.hcallback, (1, :))) + 
                soil.tillage_frac * soil.litc .* reshape(crop_cal.scallback, (1, :)) +
                (soil.tillage_frac * (soil.litc .+ max.(soil.c_input, 0.0f0))) .* reshape(crop_cal.hcallback, (1, :)) 
end

