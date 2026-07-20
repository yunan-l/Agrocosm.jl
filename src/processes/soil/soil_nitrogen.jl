"""
soil_nitrogen!(crop_cal, crop, soil)

Update litter and soil nitrogen pools and crop-available mineral nitrogen.
"""
function soil_nitrogen!(crop_cal::CropCalendar,
                        soil::Soil;
                        air_temperature = nothing,
                        wind_speed = nothing,
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack k_soil10 = lpjmlparams
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # compute soil carbon: litter carbon and soil carbon
    # soil.nitrogen.decomposed_litter = (1.0f0 .- exp.(-soil.nitrogen.litter_response / 100)) .* soil.nitrogen.litter

    soil.nitrogen.decomposed_litter = (1.0f0 .- exp.(-soil.nitrogen.litter_response .* soil.decomposition.litter_response)) .* soil.nitrogen.litter
    soil.nitrogen.litter = soil.nitrogen.litter  - soil.nitrogen.decomposed_litter

    route_harvest_nitrogen_input!(soil, crop_cal)

    litter_to_fast = soil.nitrogen.shift_fast .* sum(soil.nitrogen.decomposed_litter, dims = 1)
    litter_to_slow = soil.nitrogen.shift_slow .* sum(soil.nitrogen.decomposed_litter, dims = 1)

    # soil.nitrogen.decomposed_fast = (1.0f0 .- exp.(-soil.response_fastn .* response / 50)) .* soil.nitrogen.fast
    soil.nitrogen.decomposed_fast = (1.0f0 .- exp.(-k_soil10.fast .* soil.decomposition.response)) .* soil.nitrogen.fast
    soil.nitrogen.fast = soil.nitrogen.fast + litter_to_fast - soil.nitrogen.decomposed_fast

    # soil.nitrogen.decomposed_slow = (1.0f0 .- exp.(-soil.response_slown .* response / 10)) .* soil.nitrogen.slow
    soil.nitrogen.decomposed_slow = (1.0f0 .- exp.(-k_soil10.slow .* soil.decomposition.response)) .* soil.nitrogen.slow
    soil.nitrogen.slow = soil.nitrogen.slow + litter_to_slow - soil.nitrogen.decomposed_slow

    nitrogen_transform!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )

end
