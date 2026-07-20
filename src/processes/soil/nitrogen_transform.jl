"""
nitrogen_transform!(soil; air_temperature=nothing, wind_speed=nothing,
                    lpjmlparams=lpjmlparams)

Apply the LPJmL-style daily mineralization/immobilization, nitrification,
denitrification, and NH₃-volatilization sequence. Internal and boundary
fluxes are stored layer-wise in `soil.nitrogen` for diagnostics.
"""
function nitrogen_transform!(soil::Soil;
                             air_temperature = nothing,
                             wind_speed = nothing,
                             lpjmlparams::LPJmLParams = lpjmlparams)
    soil_layers = size(soil.nitrogen.nitrate, 1)
    decomposed_litter_carbon = vec(sum(soil.carbon.decomposed_litter; dims = 1))
    decomposed_litter_nitrogen = vec(sum(soil.nitrogen.decomposed_litter; dims = 1))

    launch_1D!(
        mineralize_immobilize_kernel!,
        decomposed_litter_carbon,
        decomposed_litter_nitrogen,
        soil.nitrogen.decomposed_fast,
        soil.nitrogen.decomposed_slow,
        soil.nitrogen.shift_fast,
        soil.nitrogen.shift_slow,
        soil.nitrogen.ammonium,
        soil.nitrogen.nitrate,
        soil.nitrogen.fast,
        soil.nitrogen.slow,
        soil.properties.layer_depth,
        soil.nitrogen.mineralization,
        soil.nitrogen.immobilization,
        (; lpjmlparams, soil_layers),
    )

    launch_1D!(
        nitrify_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        soil.nitrogen.nitrate,
        soil.water.relative_content,
        soil.water.holding_capacity_storage,
        soil.water.wilting_storage,
        soil.water.wilting_ice_fraction,
        soil.water.free_water,
        soil.water.saturation_storage,
        soil.thermal.temperature,
        soil.nitrogen.nitrification,
        soil.nitrogen.n2o_nitrification,
        (; lpjmlparams, soil_layers),
    )

    launch_1D!(
        denitrify_kernel!,
        soil.properties.ph,
        soil.carbon.fast,
        soil.carbon.slow,
        soil.water.relative_content,
        soil.water.holding_capacity_storage,
        soil.water.wilting_storage,
        soil.water.wilting_ice_fraction,
        soil.water.free_water,
        soil.water.saturation_storage,
        soil.thermal.temperature,
        soil.nitrogen.nitrate,
        soil.nitrogen.denitrification,
        soil.nitrogen.n2o_denitrification,
        soil.nitrogen.n2_denitrification,
        (; lpjmlparams, soil_layers),
    )

    volatilization_temperature = air_temperature === nothing ?
        vec(@view(soil.thermal.temperature[1, :])) : air_temperature
    volatilization_wind = if wind_speed === nothing
        fallback = similar(volatilization_temperature)
        fill!(fallback, eltype(fallback)(lpjmlparams.volatil_wind))
        fallback
    else
        wind_speed
    end
    launch_1D!(
        volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        volatilization_temperature,
        volatilization_wind,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function mineralize_immobilize_kernel!(
    decomposed_litter_carbon::AbstractArray{T},
    decomposed_litter_nitrogen::AbstractArray{T},
    decomposed_fast_nitrogen::AbstractArray{M},
    decomposed_slow_nitrogen::AbstractArray{M},
    shift_fast::AbstractArray{M},
    shift_slow::AbstractArray{M},
    ammonium::AbstractArray{M},
    nitrate::AbstractArray{M},
    fast_nitrogen::AbstractArray{M},
    slow_nitrogen::AbstractArray{M},
    layer_depth::AbstractArray{T},
    mineralization::AbstractArray{M},
    immobilization::AbstractArray{M},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack atmfrac, fastfrac, soil_cn_ratio, immobilization_k = lpjmlparams

    carbon_nitrogen_deficit =
        decomposed_litter_carbon[cell] / T(soil_cn_ratio) -
        decomposed_litter_nitrogen[cell]

    for layer in 1:soil_layers
        mineralization[layer, cell] = zero(M)
        immobilization[layer, cell] = zero(M)

        # LPJmL keeps c_shift as a normalized vertical distribution. Apply the
        # fast/slow split and atmospheric fraction explicitly to each flux.
        litter_mineralization = max(
            zero(M),
            decomposed_litter_nitrogen[cell] * T(atmfrac) *
            (T(fastfrac) * shift_fast[layer, cell] +
             (one(T) - T(fastfrac)) * shift_slow[layer, cell]),
        )
        som_mineralization = max(
            zero(M),
            decomposed_fast_nitrogen[layer, cell] +
            decomposed_slow_nitrogen[layer, cell],
        )
        gross_mineralization = litter_mineralization + som_mineralization
        ammonium[layer, cell] += gross_mineralization
        mineralization[layer, cell] = gross_mineralization

        if carbon_nitrogen_deficit > zero(T)
            available = max(zero(M), ammonium[layer, cell] + nitrate[layer, cell])
            if available > zero(M)
                concentration = available / max(layer_depth[layer], eps(T)) * T(1000)
                limitation = concentration / (T(immobilization_k) + concentration)

                fast_immobilization = max(
                    zero(M),
                    carbon_nitrogen_deficit * T(fastfrac) *
                    (one(T) - T(atmfrac)) * shift_fast[layer, cell] * limitation,
                )
                fast_immobilization = min(fast_immobilization, available)
                if fast_immobilization > zero(M)
                    ammonium_share = ammonium[layer, cell] / available
                    ammonium[layer, cell] -= fast_immobilization * ammonium_share
                    nitrate[layer, cell] -= fast_immobilization * (one(M) - ammonium_share)
                    fast_nitrogen[layer, cell] += fast_immobilization
                    immobilization[layer, cell] += fast_immobilization
                end

                available = max(zero(M), ammonium[layer, cell] + nitrate[layer, cell])
                if available > zero(M)
                    concentration = available / max(layer_depth[layer], eps(T)) * T(1000)
                    limitation = concentration / (T(immobilization_k) + concentration)
                    slow_immobilization = max(
                        zero(M),
                        carbon_nitrogen_deficit * (one(T) - T(fastfrac)) *
                        (one(T) - T(atmfrac)) * shift_slow[layer, cell] * limitation,
                    )
                    slow_immobilization = min(slow_immobilization, available)
                    if slow_immobilization > zero(M)
                        ammonium_share = ammonium[layer, cell] / available
                        ammonium[layer, cell] -= slow_immobilization * ammonium_share
                        nitrate[layer, cell] -= slow_immobilization * (one(M) - ammonium_share)
                        slow_nitrogen[layer, cell] += slow_immobilization
                        immobilization[layer, cell] += slow_immobilization
                    end
                end
            end
        end
        ammonium[layer, cell] = max(zero(M), ammonium[layer, cell])
        nitrate[layer, cell] = max(zero(M), nitrate[layer, cell])
    end
end

@kernel inbounds = true function nitrify_kernel!(
    soil_ph::AbstractArray{T},
    ammonium::AbstractArray{M},
    nitrate::AbstractArray{M},
    relative_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    wilting_storage::AbstractArray{M},
    wilting_ice_fraction::AbstractArray{M},
    free_water::AbstractArray{M},
    saturation_storage::AbstractArray{M},
    soil_temperature::AbstractArray{M},
    nitrification::AbstractArray{M},
    n2o_nitrification::AbstractArray{M},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack k_max, k_2, nitrification_a, nitrification_b,
            nitrification_c, nitrification_d = lpjmlparams

    for layer in 1:soil_layers
        nitrification[layer, cell] = zero(M)
        n2o_nitrification[layer, cell] = zero(M)
        water_filled_pore_space = clamp(
            (relative_water[layer, cell] * holding_storage[layer, cell] +
             wilting_storage[layer, cell] *
             (one(M) - wilting_ice_fraction[layer, cell]) +
             free_water[layer, cell]) /
            max(saturation_storage[layer, cell], eps(M)),
            zero(M), one(M),
        )
        n_nit = T(nitrification_a - nitrification_b)
        m_nit = T(nitrification_a - nitrification_c)
        z_nit = T(nitrification_d) * T(nitrification_b - nitrification_a) /
            T(nitrification_a - nitrification_c)
        base_1 = (water_filled_pore_space - T(nitrification_b)) / n_nit
        base_2 = (water_filled_pore_space - T(nitrification_c)) / m_nit
        moisture_factor = if base_1 > zero(M) && base_2 > zero(M)
            max(zero(M), base_1^z_nit * base_2^T(nitrification_d))
        else
            zero(M)
        end
        temperature_factor = exp(
            -(soil_temperature[layer, cell] - T(18.79))^2 /
            T(2 * 8.26 * 8.26),
        )
        ph_factor = T(0.56) +
            atan(T(pi) * T(0.45) * (soil_ph[cell] - T(5))) / T(pi)
        gross_nitrification = clamp(
            T(k_max) * ammonium[layer, cell] * temperature_factor *
            moisture_factor * ph_factor,
            zero(M), ammonium[layer, cell],
        )
        n2o_loss = T(k_2) * gross_nitrification
        ammonium[layer, cell] -= gross_nitrification
        nitrate[layer, cell] += gross_nitrification - n2o_loss
        nitrification[layer, cell] = gross_nitrification
        n2o_nitrification[layer, cell] = n2o_loss
    end
end

@kernel inbounds = true function denitrify_kernel!(
    _soil_ph::AbstractArray{T},
    fast_carbon::AbstractArray{M},
    slow_carbon::AbstractArray{M},
    relative_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    wilting_storage::AbstractArray{M},
    wilting_ice_fraction::AbstractArray{M},
    free_water::AbstractArray{M},
    saturation_storage::AbstractArray{M},
    soil_temperature::AbstractArray{M},
    nitrate::AbstractArray{M},
    denitrification::AbstractArray{M},
    n2o_denitrification::AbstractArray{M},
    n2_denitrification::AbstractArray{M},
    kernel_params,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack CDN, n2o_denit_frac = lpjmlparams

    for layer in 1:soil_layers
        denitrification[layer, cell] = zero(M)
        n2o_denitrification[layer, cell] = zero(M)
        n2_denitrification[layer, cell] = zero(M)
        temperature = soil_temperature[layer, cell]
        organic_carbon = max(zero(M), fast_carbon[layer, cell] + slow_carbon[layer, cell])
        temperature_factor = if temperature > M(45.9)
            zero(M)
        elseif temperature > zero(M)
            max(zero(M), M(0.0326) + M(0.00351) * temperature^M(1.652) -
                (temperature / M(41.748))^M(7.19))
        else
            M(0.0326)
        end
        water_filled_pore_space =
            (wilting_storage[layer, cell] *
             (one(M) - wilting_ice_fraction[layer, cell]) +
             relative_water[layer, cell] * holding_storage[layer, cell] +
             free_water[layer, cell]) /
            max(saturation_storage[layer, cell], eps(M))
        moisture_factor = min(
            one(M), M(6.664096e-10) * exp(M(20.92912) * water_filled_pore_space),
        )
        carbon_factor = max(
            zero(M), one(M) - exp(-M(CDN) * temperature_factor * organic_carbon),
        )
        gross_denitrification = clamp(
            moisture_factor * carbon_factor * nitrate[layer, cell],
            zero(M), nitrate[layer, cell],
        )
        n2o_loss = M(n2o_denit_frac) * gross_denitrification
        n2_loss = gross_denitrification - n2o_loss
        nitrate[layer, cell] -= gross_denitrification
        denitrification[layer, cell] = gross_denitrification
        n2o_denitrification[layer, cell] = n2o_loss
        n2_denitrification[layer, cell] = n2_loss
    end
end

@kernel inbounds = true function volatilization_kernel!(
    soil_ph::AbstractArray{T},
    ammonium::AbstractArray{M},
    air_temperature::AbstractArray{M},
    wind_speed::AbstractArray{M},
    layer_depth::AbstractArray{T},
    volatilization::AbstractArray{M},
    lpjmlparams::LPJmLParams,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    @unpack volatil_length = lpjmlparams
    temperature = air_temperature[cell]
    kelvin = temperature + M(273.15)
    ammonium_top = max(zero(M), ammonium[1, cell])
    dissociation = M(10)^(M(0.05) - M(2788) / kelvin)
    aqueous_fraction = one(M) /
        (one(M) + M(10)^(-soil_ph[cell]) / max(dissociation, eps(M)))
    aqueous_nh3 = aqueous_fraction * ammonium_top /
        max(layer_depth[1], eps(T)) * M(1000)
    henry = M(0.2138) / kelvin * M(10)^(M(6.123) - M(1825) / kelvin)
    mass_transfer = M(0.000612) * max(zero(M), wind_speed[cell])^M(0.8) *
        kelvin^M(0.382) * M(volatil_length)^M(-0.2)
    flux = clamp(M(86400) * mass_transfer * henry * aqueous_nh3,
                 zero(M), ammonium_top)
    ammonium[1, cell] -= flux
    volatilization[cell] = flux
end
