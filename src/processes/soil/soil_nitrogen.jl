"""
soil_nitrogen!(crop_cal, crop, soil)

Update litter and soil nitrogen pools and crop-available mineral nitrogen.
"""
function soil_nitrogen!(crop_cal::CropCalendar,
                        soil::Soil;
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack k_soil10 = lpjmlparams
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # compute soil carbon: litter carbon and soil carbon
    # soil.nitrogen.decomposed_litter = (1.0f0 .- exp.(-soil.nitrogen.litter_response / 100)) .* soil.nitrogen.litter

    soil.nitrogen.decomposed_litter = (1.0f0 .- exp.(-soil.nitrogen.litter_response .* soil.decomposition.litter_response)) .* soil.nitrogen.litter
    soil.nitrogen.litter = soil.nitrogen.litter  - soil.nitrogen.decomposed_litter

    # using 'callback' to adjust litter carbon due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    update_litn_tillage!(soil, crop_cal)

    litter_to_fast = soil.nitrogen.shift_fast .* sum(soil.nitrogen.decomposed_litter, dims = 1)
    litter_to_slow = soil.nitrogen.shift_slow .* sum(soil.nitrogen.decomposed_litter, dims = 1)

    # soil.nitrogen.decomposed_fast = (1.0f0 .- exp.(-soil.response_fastn .* response / 50)) .* soil.nitrogen.fast
    soil.nitrogen.decomposed_fast = (1.0f0 .- exp.(-k_soil10.fast .* soil.decomposition.response)) .* soil.nitrogen.fast
    soil.nitrogen.fast = soil.nitrogen.fast + litter_to_fast - soil.nitrogen.decomposed_fast

    # soil.nitrogen.decomposed_slow = (1.0f0 .- exp.(-soil.response_slown .* response / 10)) .* soil.nitrogen.slow
    soil.nitrogen.decomposed_slow = (1.0f0 .- exp.(-k_soil10.slow .* soil.decomposition.response)) .* soil.nitrogen.slow
    soil.nitrogen.slow = soil.nitrogen.slow + litter_to_slow - soil.nitrogen.decomposed_slow

    # Conservative M1 mineralization baseline. Litter N not retained in the
    # fast/slow organic pools, plus decomposed fast/slow N, becomes NH4. The
    # M3 transformation module will later partition this mineral N and record
    # nitrification, denitrification, volatilization, and leaching explicitly.
    mineralized_litter = max.(
        zero(eltype(soil.nitrogen.ammonium)),
        vec(sum(soil.nitrogen.decomposed_litter; dims = 1)) .-
        vec(sum(litter_to_fast .+ litter_to_slow; dims = 1)),
    )
    @views soil.nitrogen.ammonium[1, :] .+= mineralized_litter
    soil.nitrogen.ammonium .+=
        soil.nitrogen.decomposed_fast .+ soil.nitrogen.decomposed_slow

    # Nitrogen mineralization/immobilization/nitrification updates.
    # TODO: we do not yet test it, so we close it for now.
    # nitrogen_transform!(soil; lpjmlparams = lpjmlparams)

end


"""
update_litn_tillage!(soil, crop_cal)

Apply tillage/harvest crop nitrogen to litter nitrogen pools.
"""
function update_litn_tillage!(soil::Soil,
                              crop_cal::CropCalendar
)

    soil.nitrogen.litter = soil.nitrogen.litter .* (1 .- reshape(crop_cal.sowing_callback, (1, :))) .* (1 .- reshape(crop_cal.harvest_callback, (1, :))) +
                soil.management.tillage_fraction * soil.nitrogen.litter .* reshape(crop_cal.sowing_callback, (1, :)) +
                (soil.management.tillage_fraction * (soil.nitrogen.litter .+ max.(soil.nitrogen.input, 0.0f0))) .* reshape(crop_cal.harvest_callback, (1, :))

end
