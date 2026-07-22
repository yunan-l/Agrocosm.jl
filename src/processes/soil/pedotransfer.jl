"""
pedotransfer!(soil; lpjmlparams=lpjmlparams)

Derive soil hydraulic properties from texture and depth parameterizations.
"""
function pedotransfer_reference!(soil;
                                 lpjmlparams::LPJmLParams = lpjmlparams
)

    @unpack MINERALDENS = lpjmlparams # mineral density in kg/m3

    om_layer = 2 * ((soil_carbon_prognostic(soil).fast + soil_carbon_prognostic(soil).slow) ./ ((1 .- soil_water_prognostic(soil).saturation_fraction) * MINERALDENS .* soil_properties(soil).layer_depth)) * 100 #calculation of soil organic matter in %

    # idx = om_layer .> 8
    # om_layer[idx] .= T(8.0)
    # om_layer .= ifelse.(om_layer .> 8, 8.0, om_layer)
    om_layer .= clamp.(om_layer, 0.0f0, 8.0f0)

    wpwpt = -0.024f0 * soil_properties(soil).sand_fraction + 0.487f0 * soil_properties(soil).clay_fraction .+ 0.006f0 * om_layer + 0.005f0 * (soil_properties(soil).sand_fraction .* om_layer) - 0.013f0 * (soil_properties(soil).clay_fraction .* om_layer) .+ 0.068f0 * (soil_properties(soil).sand_fraction .* soil_properties(soil).clay_fraction) .+ 0.031f0
    soil_water_auxiliary(soil).wilting_fraction .= wpwpt + (0.14 * wpwpt .- 0.02)
    soil_water_auxiliary(soil).wilting_storage .= soil_water_auxiliary(soil).wilting_fraction .* soil_properties(soil).layer_depth
    ws33t = 0.278f0 * soil_properties(soil).sand_fraction + 0.034f0 * soil_properties(soil).clay_fraction .+ 0.022f0 * om_layer - 0.018f0 * (soil_properties(soil).sand_fraction .* om_layer) - 0.027f0 * (soil_properties(soil).clay_fraction .* om_layer) .- 0.584f0 * (soil_properties(soil).sand_fraction .* soil_properties(soil).clay_fraction) .+ 0.078f0
    ws33 = ws33t + (0.636f0 * ws33t .- 0.107f0)

    wfct = -0.251f0 * soil_properties(soil).sand_fraction + 0.195f0 * soil_properties(soil).clay_fraction .+ 0.011f0 * om_layer + 0.006f0 * (soil_properties(soil).sand_fraction .* om_layer) - 0.027f0 * (soil_properties(soil).clay_fraction .* om_layer) .+ 0.452f0 * (soil_properties(soil).sand_fraction .* soil_properties(soil).clay_fraction) .+ 0.299f0
    soil_water_auxiliary(soil).field_capacity .= (wfct + (((1.283f0 * wfct) .^ 2) - 0.374f0 * wfct .- 0.015f0))

    base_saturation = soil_water_auxiliary(soil).field_capacity + ws33 .-
        0.097f0 * soil_properties(soil).sand_fraction .+ 0.043f0
    soil_water_prognostic(soil).saturation_fraction .= base_saturation

    # LPJmL tills only the upper soil layer. Lower bulk density increases its
    # pore volume and slightly shifts field capacity until rainfall settles it.
    @views begin
        density_factor = soil_management_prognostic(soil).tillage_density_factor[1, :]
        tilled_saturation = one(eltype(base_saturation)) .-
            (one(eltype(base_saturation)) .- base_saturation[1, :]) .* density_factor
        soil_water_prognostic(soil).saturation_fraction[1, :] .= tilled_saturation
        soil_water_auxiliary(soil).field_capacity[1, :] .-= 0.2f0 .* (
            base_saturation[1, :] .- tilled_saturation
        )
    end
    soil_water_auxiliary(soil).saturation_storage .= soil_water_prognostic(soil).saturation_fraction .* soil_properties(soil).layer_depth

    # idx = (soil_water_prognostic(soil).saturation_fraction - soil_water_auxiliary(soil).field_capacity) .< 0.05
    # soil_water_auxiliary(soil).field_capacity[idx] .= soil_water_prognostic(soil).saturation_fraction[idx] .- 0.05
    soil_water_auxiliary(soil).field_capacity .= ifelse.((soil_water_prognostic(soil).saturation_fraction - soil_water_auxiliary(soil).field_capacity) .< 0.05f0, soil_water_prognostic(soil).saturation_fraction .- 0.05f0, soil_water_auxiliary(soil).field_capacity)

    soil_water_auxiliary(soil).beta .= -2.655f0 ./ log10.(soil_water_auxiliary(soil).field_capacity ./ soil_water_prognostic(soil).saturation_fraction)
    soil_water_auxiliary(soil).holding_capacity_fraction .= soil_water_auxiliary(soil).field_capacity - soil_water_auxiliary(soil).wilting_fraction
    soil_water_auxiliary(soil).holding_capacity_storage .= soil_water_auxiliary(soil).holding_capacity_fraction .* soil_properties(soil).layer_depth

    # Reconstruct LPJmL's below-PWP, plant-available, and free-water pools
    # after hydraulic capacities change. Total liquid water and total ice are
    # the conserved states; the individual reservoirs are deterministic.
    partition_soil_water_ice!(soil)

    # Calculation of Ks
    lambda = (log.(soil_water_auxiliary(soil).field_capacity) - log.(soil_water_auxiliary(soil).wilting_fraction)) / (log(1500) - log(33))
    soil_water_auxiliary(soil).saturated_conductivity .= 1930 * (soil_water_prognostic(soil).saturation_fraction - soil_water_auxiliary(soil).field_capacity) .^ (3 .- lambda)

end

"""Compute hydraulic properties in one backend-neutral layer/cell kernel."""
function pedotransfer!(soil;
                       lpjmlparams::LPJmLParams = lpjmlparams)
    launch_2D!(
        pedotransfer_kernel!,
        soil_water_auxiliary(soil).wilting_fraction,
        soil_properties(soil).sand_fraction,
        soil_properties(soil).clay_fraction,
        soil_properties(soil).layer_depth,
        soil_carbon_prognostic(soil).fast,
        soil_carbon_prognostic(soil).slow,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_auxiliary(soil).field_capacity,
        soil_water_prognostic(soil).saturation_fraction,
        soil_water_auxiliary(soil).saturation_storage,
        soil_water_auxiliary(soil).beta,
        soil_water_auxiliary(soil).holding_capacity_fraction,
        soil_water_auxiliary(soil).holding_capacity_storage,
        soil_water_auxiliary(soil).saturated_conductivity,
        soil_management_prognostic(soil).tillage_density_factor,
        eltype(soil_water_prognostic(soil).storage)(lpjmlparams.MINERALDENS),
    )
    partition_soil_water_ice!(soil)
    return nothing
end

function pedotransfer!(state::ModelState;
                       lpjmlparams::LPJmLParams = lpjmlparams)
    water_state = state.prognostic.soil.water
    water_auxiliary = state.auxiliary.soil.water
    soil_inputs = state.inputs.soil
    soil_carbon = state.prognostic.soil.carbon
    launch_2D!(
        pedotransfer_kernel!, water_auxiliary.wilting_fraction,
        soil_inputs.properties.sand_fraction,
        soil_inputs.properties.clay_fraction,
        soil_inputs.properties.layer_depth, soil_carbon.fast, soil_carbon.slow,
        water_auxiliary.wilting_storage, water_auxiliary.field_capacity,
        water_state.saturation_fraction, water_auxiliary.saturation_storage,
        water_auxiliary.beta, water_auxiliary.holding_capacity_fraction,
        water_auxiliary.holding_capacity_storage,
        water_auxiliary.saturated_conductivity,
        state.prognostic.soil.management.tillage_density_factor,
        eltype(water_state.storage)(lpjmlparams.MINERALDENS),
    )
    partition_soil_water_ice!(state)
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
    tillage_density_factor::AbstractMatrix{T},
    mineral_density::T,
) where {T <: AbstractFloat}
    layer, cell = @index(Global, NTuple)

    sand = sand_fraction[1, cell]
    clay = clay_fraction[1, cell]
    depth = layer_depth[layer]
    previous_saturation = saturation_fraction[layer, cell]
    organic_matter = clamp(
        T(2) * ((fast_carbon[layer, cell] + slow_carbon[layer, cell]) /
        ((one(T) - previous_saturation) * mineral_density * depth)) * T(100),
        zero(T),
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
    base_saturation = field + ws33 - T(0.097) * sand + T(0.043)
    if layer == 1
        saturation = one(T) -
            (one(T) - base_saturation) * tillage_density_factor[1, cell]
        field -= T(0.2) * (base_saturation - saturation)
    else
        saturation = base_saturation
    end
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
