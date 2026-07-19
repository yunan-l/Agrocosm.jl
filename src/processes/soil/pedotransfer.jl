"""
pedotransfer!(soil; lpjmlparams=lpjmlparams)

Derive soil hydraulic properties from texture and depth parameterizations.
"""
function pedotransfer!(soil::Soil;
                       lpjmlparams::LPJmLParams = lpjmlparams
)

    @unpack MINERALDENS = lpjmlparams # mineral density in kg/m3

    om_layer = 2 * ((soil.carbon.fast + soil.carbon.slow) ./ ((1 .- soil.water.saturation_fraction) * MINERALDENS .* soil.properties.layer_depth)) * 100 #calculation of soil organic matter in %

    # idx = om_layer .> 8
    # om_layer[idx] .= T(8.0)
    # om_layer .= ifelse.(om_layer .> 8, 8.0, om_layer)
    om_layer .= min.(om_layer, 8.0f0)

    wpwpt = -0.024f0 * soil.properties.sand_fraction + 0.487f0 * soil.properties.clay_fraction .+ 0.006f0 * om_layer + 0.005f0 * (soil.properties.sand_fraction .* om_layer) - 0.013f0 * (soil.properties.clay_fraction .* om_layer) .+ 0.068f0 * (soil.properties.sand_fraction .* soil.properties.clay_fraction) .+ 0.031f0
    soil.water.wilting_fraction .= wpwpt + (0.14 * wpwpt .- 0.02)
    soil.water.wilting_storage .= soil.water.wilting_fraction .* soil.properties.layer_depth
    ws33t = 0.278f0 * soil.properties.sand_fraction + 0.034f0 * soil.properties.clay_fraction .+ 0.022f0 * om_layer - 0.018f0 * (soil.properties.sand_fraction .* om_layer) - 0.027f0 * (soil.properties.clay_fraction .* om_layer) .- 0.584f0 * (soil.properties.sand_fraction .* soil.properties.clay_fraction) .+ 0.078f0
    ws33 = ws33t + (0.636f0 * ws33t .- 0.107f0)

    wfct = -0.251f0 * soil.properties.sand_fraction + 0.195f0 * soil.properties.clay_fraction .+ 0.011f0 * om_layer + 0.006f0 * (soil.properties.sand_fraction .* om_layer) - 0.027f0 * (soil.properties.clay_fraction .* om_layer) .+ 0.452f0 * (soil.properties.sand_fraction .* soil.properties.clay_fraction) .+ 0.299f0
    soil.water.field_capacity .= (wfct + (((1.283f0 * wfct) .^ 2) - 0.374f0 * wfct .- 0.015f0))

    soil.water.saturation_fraction .= soil.water.field_capacity + ws33 .- 0.097f0 * soil.properties.sand_fraction .+ 0.043f0
    soil.water.saturation_storage .= soil.water.saturation_fraction .* soil.properties.layer_depth

    # here, we ignore the effects of tillage to soil water content at saturation.
    # if(l < NTILLLAYER)
    # {
    #     soil->wsat[l] = 1 - (1-w_sat)*soil->df_tillage[l];
    #     soil->wfc[l] = w_fc - 0.2 * (w_sat - soil->wsat[l]);
    # }
    # else
    # {
    #     soil->wsat[l] = w_sat;
    #     soil->wfc[l] = w_fc;
    # }

    # idx = (soil.water.saturation_fraction - soil.water.field_capacity) .< 0.05
    # soil.water.field_capacity[idx] .= soil.water.saturation_fraction[idx] .- 0.05
    soil.water.field_capacity .= ifelse.((soil.water.saturation_fraction - soil.water.field_capacity) .< 0.05f0, soil.water.saturation_fraction .- 0.05f0, soil.water.field_capacity)

    soil.water.beta .= -2.655f0 ./ log10.(soil.water.field_capacity ./ soil.water.saturation_fraction)
    soil.water.holding_capacity_fraction .= soil.water.field_capacity - soil.water.wilting_fraction
    soil.water.holding_capacity_storage .= soil.water.holding_capacity_fraction .* soil.properties.layer_depth

    # idx = (soil.water.storage - soil.water.wilting_storage) .> 1.0f-10
    # soil.water.relative_content[idx] .= min.((soil.water.storage[idx] - soil.water.wilting_storage[idx]) ./ soil.water.holding_capacity_storage[idx], one(T))
    # soil.water.free_water[idx] .= min.(soil.water.storage[idx] - soil.water.wilting_storage[idx] - soil.water.relative_content[idx] .* soil.water.holding_capacity_storage[idx], soil.water.saturation_storage[idx] - soil.water.field_capacity[idx] .* soil.properties.layer_depth)
    # idx = (soil.water.storage - soil.water.wilting_storage) <= 1.0f-10
    # soil.water.relative_content[idx] .= zero(T)
    # soil.water.free_water[idx] .= zero(T)
    soil.water.relative_content .= ifelse.((soil.water.storage - soil.water.wilting_storage) .> 1.0f-10, min.((soil.water.storage - soil.water.wilting_storage) ./ soil.water.holding_capacity_storage, 1.0f0), 0.0f0)
    soil.water.free_water .= ifelse.((soil.water.storage - soil.water.wilting_storage) .> 1.0f-10, min.(soil.water.storage - soil.water.wilting_storage - soil.water.relative_content .* soil.water.holding_capacity_storage, soil.water.saturation_storage - soil.water.field_capacity .* soil.properties.layer_depth), 0.0f0)

    # Calculation of Ks
    lambda = (log.(soil.water.field_capacity) - log.(soil.water.wilting_fraction)) / (log(1500) - log(33))
    soil.water.saturated_conductivity .= 1930 * (soil.water.saturation_fraction - soil.water.field_capacity) .^ (3 .- lambda)

    # update agtop_cover
    dm_sum = soil.carbon.litter[1, :]/0.42 # dm_sum+=stand->soil.litter.item[l].agtop.leaf.carbon/0.42; Accounting that C content in plant dry matter is 42%
    # idx = dm_sum .< 0
    # dm_sum[idx] .= zero(T)
    # dm_sum .= ifelse.(dm_sum .< 0, 0.0, dm_sum)
    dm_sum .= max.(dm_sum, 0.0f0)
    soil.properties.surface_litter_cover .= 1 .- exp.(-6e-3 * dm_sum)
end
