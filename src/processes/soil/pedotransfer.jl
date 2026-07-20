"""
pedotransfer!(soil; lpjmlparams=lpjmlparams)

Derive soil hydraulic properties from texture and depth parameterizations.
"""
function pedotransfer_reference!(soil::Soil;
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

    # Reconstruct LPJmL's below-PWP, plant-available, and free-water pools
    # after hydraulic capacities change. Total liquid water and total ice are
    # the conserved states; the individual reservoirs are deterministic.
    partition_soil_water_ice!(soil)

    # Calculation of Ks
    lambda = (log.(soil.water.field_capacity) - log.(soil.water.wilting_fraction)) / (log(1500) - log(33))
    soil.water.saturated_conductivity .= 1930 * (soil.water.saturation_fraction - soil.water.field_capacity) .^ (3 .- lambda)

end

"""Compute hydraulic properties in one backend-neutral layer/cell kernel."""
function pedotransfer!(soil::Soil;
                       lpjmlparams::LPJmLParams = lpjmlparams)
    launch_2D!(
        pedotransfer_kernel!,
        soil.water.wilting_fraction,
        soil.properties.sand_fraction,
        soil.properties.clay_fraction,
        soil.properties.layer_depth,
        soil.carbon.fast,
        soil.carbon.slow,
        soil.water.wilting_storage,
        soil.water.field_capacity,
        soil.water.saturation_fraction,
        soil.water.saturation_storage,
        soil.water.beta,
        soil.water.holding_capacity_fraction,
        soil.water.holding_capacity_storage,
        soil.water.saturated_conductivity,
        eltype(soil.water.storage)(lpjmlparams.MINERALDENS),
    )
    partition_soil_water_ice!(soil)
    return nothing
end

@kernel inbounds = true function pedotransfer_kernel!(
    wilting_fraction::AbstractMatrix{T},
    sand_fraction::AbstractMatrix{T},
    clay_fraction::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    fast_carbon::AbstractMatrix{T},
    slow_carbon::AbstractMatrix{T},
    wilting_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    saturation_fraction::AbstractMatrix{T},
    saturation_storage::AbstractMatrix{T},
    beta::AbstractMatrix{T},
    holding_capacity_fraction::AbstractMatrix{T},
    holding_capacity_storage::AbstractMatrix{T},
    saturated_conductivity::AbstractMatrix{T},
    mineral_density::T,
) where {T <: AbstractFloat}
    layer, cell = @index(Global, NTuple)

    sand = sand_fraction[1, cell]
    clay = clay_fraction[1, cell]
    depth = layer_depth[layer]
    previous_saturation = saturation_fraction[layer, cell]
    organic_matter = min(
        T(2) * ((fast_carbon[layer, cell] + slow_carbon[layer, cell]) /
        ((one(T) - previous_saturation) * mineral_density * depth)) * T(100),
        T(8),
    )

    wpwpt = -T(0.024) * sand + T(0.487) * clay + T(0.006) * organic_matter +
        T(0.005) * sand * organic_matter - T(0.013) * clay * organic_matter +
        T(0.068) * sand * clay + T(0.031)
    wilting = wpwpt + (T(0.14) * wpwpt - T(0.02))

    ws33t = T(0.278) * sand + T(0.034) * clay + T(0.022) * organic_matter -
        T(0.018) * sand * organic_matter - T(0.027) * clay * organic_matter -
        T(0.584) * sand * clay + T(0.078)
    ws33 = ws33t + (T(0.636) * ws33t - T(0.107))

    wfct = -T(0.251) * sand + T(0.195) * clay + T(0.011) * organic_matter +
        T(0.006) * sand * organic_matter - T(0.027) * clay * organic_matter +
        T(0.452) * sand * clay + T(0.299)
    field = wfct + ((T(1.283) * wfct)^2 - T(0.374) * wfct - T(0.015))
    saturation = field + ws33 - T(0.097) * sand + T(0.043)
    field = saturation - field < T(0.05) ? saturation - T(0.05) : field

    wilting_fraction[layer, cell] = wilting
    wilting_storage[layer, cell] = wilting * depth
    field_capacity[layer, cell] = field
    saturation_fraction[layer, cell] = saturation
    saturation_storage[layer, cell] = saturation * depth
    beta[layer, cell] = -T(2.655) / log10(field / saturation)
    holding = field - wilting
    holding_capacity_fraction[layer, cell] = holding
    holding_capacity_storage[layer, cell] = holding * depth
    lambda = (log(field) - log(wilting)) / (log(T(1500)) - log(T(33)))
    saturated_conductivity[layer, cell] =
        T(1930) * (saturation - field)^(T(3) - lambda)
end
