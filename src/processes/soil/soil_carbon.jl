"""
soil_carbon!(crop_cal, crop, soil)

Update litter and soil carbon pools and heterotrophic respiration terms.
"""
function soil_carbon!(crop_cal::CropCalendar,
                      soil::Soil;
                      lpjmlparams::LPJmLParams = lpjmlparams,
                      soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack atmfrac, fastfrac, k_soil10 = lpjmlparams
    @unpack e0, intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # soil decomposition response
    soil_decomp_response!(soil)

    # compute soil carbon: litter carbon and soil carbon
    # soil.carbon.decomposed_litter = (1.0f0 .- exp.(-soil.carbon.litter_response / 100)) .* soil.carbon.litter

    # Litter decomposition is represented as three aggregated litter pools.
    # We use top-layer response (LPJmL uses top/root layer litter environments).
    soil.carbon.decomposed_litter = (1.0f0 .- exp.(-soil.carbon.litter_response .* soil.decomposition.litter_response)) .* soil.carbon.litter
    soil.carbon.litter = soil.carbon.litter  - soil.carbon.decomposed_litter

    # LPJmL harvest first creates agtop/bg litter, then the KILL -> setaside
    # transition tills agtop into agsub on the same day.
    route_harvest_carbon_input!(soil, crop_cal)

    # soil.carbon.decomposed_fast = (1.0f0 .- exp.(-soil.response_fastc .* response / 50)) .* soil.carbon.fast
    soil.carbon.decomposed_fast = (1.0f0 .- exp.(-k_soil10.fast .* soil.decomposition.response)) .* soil.carbon.fast
    soil.carbon.litter_to_fast .= soil.carbon.shift_fast .*
        sum(soil.carbon.decomposed_litter, dims = 1) .* fastfrac .* (1.0f0 - atmfrac)
    soil.carbon.fast = soil.carbon.fast + soil.carbon.litter_to_fast - soil.carbon.decomposed_fast

    # soil.carbon.decomposed_slow = (1.0f0 .- exp.(-soil.response_slowc .* response / 10)) .* soil.carbon.slow
    soil.carbon.decomposed_slow = (1.0f0 .- exp.(-k_soil10.slow .* soil.decomposition.response)) .* soil.carbon.slow
    soil.carbon.litter_to_slow .= soil.carbon.shift_slow .*
        sum(soil.carbon.decomposed_litter, dims = 1) .* (1.0f0 - fastfrac) .* (1.0f0 - atmfrac)
    soil.carbon.slow = soil.carbon.slow + soil.carbon.litter_to_slow - soil.carbon.decomposed_slow

    soil.carbon.heterotrophic_respiration = vec(sum(soil.carbon.decomposed_litter, dims = 1) * atmfrac .+ sum(soil.carbon.decomposed_fast, dims = 1) .+ sum(soil.carbon.decomposed_slow, dims = 1))

end
