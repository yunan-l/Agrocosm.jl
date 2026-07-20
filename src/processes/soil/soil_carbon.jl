"""
soil_carbon!(crop_cal, crop, soil)

Update litter and soil carbon pools and heterotrophic respiration terms.
"""
function soil_carbon_reference!(crop_cal::CropCalendar,
                                soil::Soil;
                                lpjmlparams::LPJmLParams = lpjmlparams,
                                soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack atmfrac, fastfrac, k_soil10 = lpjmlparams
    @unpack e0, intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # soil decomposition response
    soil_decomp_response_reference!(soil)

    # compute soil carbon: litter carbon and soil carbon
    # soil.carbon.decomposed_litter = (1.0f0 .- exp.(-soil.carbon.litter_response / 100)) .* soil.carbon.litter

    # Litter decomposition is represented as three aggregated litter pools.
    # We use top-layer response (LPJmL uses top/root layer litter environments).
    soil.carbon.decomposed_litter .=
        (1.0f0 .- exp.(-soil.carbon.litter_response .* soil.decomposition.litter_response)) .* soil.carbon.litter
    soil.carbon.litter .-= soil.carbon.decomposed_litter

    # LPJmL harvest first creates agtop/bg litter, then the KILL -> setaside
    # transition tills agtop into agsub on the same day.
    route_harvest_carbon_input_reference!(soil, crop_cal)

    # soil.carbon.decomposed_fast = (1.0f0 .- exp.(-soil.response_fastc .* response / 50)) .* soil.carbon.fast
    soil.carbon.decomposed_fast .= max.(
        0.0f0,
        (1.0f0 .- exp.(-k_soil10.fast .* soil.decomposition.response)) .* soil.carbon.fast,
    )
    decomposed_litter = soil.decomposition.surface_scratch_1
    @views decomposed_litter .=
        soil.carbon.decomposed_litter[1, :] .+
        soil.carbon.decomposed_litter[2, :] .+
        soil.carbon.decomposed_litter[3, :]
    soil.carbon.litter_to_fast .= soil.carbon.shift_fast .*
        reshape(decomposed_litter, 1, :) .* fastfrac .* (1.0f0 - atmfrac)
    soil.carbon.fast .+= soil.carbon.litter_to_fast .- soil.carbon.decomposed_fast

    # soil.carbon.decomposed_slow = (1.0f0 .- exp.(-soil.response_slowc .* response / 10)) .* soil.carbon.slow
    soil.carbon.decomposed_slow .= max.(
        0.0f0,
        (1.0f0 .- exp.(-k_soil10.slow .* soil.decomposition.response)) .* soil.carbon.slow,
    )
    soil.carbon.litter_to_slow .= soil.carbon.shift_slow .*
        reshape(decomposed_litter, 1, :) .* (1.0f0 - fastfrac) .* (1.0f0 - atmfrac)
    soil.carbon.slow .+= soil.carbon.litter_to_slow .- soil.carbon.decomposed_slow

    soil.carbon.heterotrophic_respiration .= decomposed_litter .* atmfrac
    for layer in axes(soil.carbon.decomposed_fast, 1)
        @views soil.carbon.heterotrophic_respiration .+=
            soil.carbon.decomposed_fast[layer, :] .+
            soil.carbon.decomposed_slow[layer, :]
    end

end

function soil_carbon!(crop_cal::CropCalendar,
                      soil::Soil;
                      lpjmlparams::LPJmLParams = lpjmlparams,
                      soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_decomp_response!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    T = eltype(soil.carbon.litter)
    launch_custom!(
        soil_carbon_decomposition_kernel!,
        soil.carbon.litter,
        size(soil.carbon.litter, 2),
        soil.carbon.litter_response,
        soil.decomposition.litter_response,
        soil.carbon.decomposed_litter,
        soil.carbon.fast,
        soil.carbon.slow,
        soil.decomposition.response,
        soil.carbon.decomposed_fast,
        soil.carbon.decomposed_slow,
        soil.carbon.shift_fast,
        soil.carbon.shift_slow,
        soil.carbon.litter_to_fast,
        soil.carbon.litter_to_slow,
        soil.carbon.heterotrophic_respiration,
        T(lpjmlparams.atmfrac),
        T(lpjmlparams.fastfrac),
        T(lpjmlparams.k_soil10.fast),
        T(lpjmlparams.k_soil10.slow),
        size(soil.carbon.fast, 1),
    )
    route_harvest_carbon_input!(soil, crop_cal)
    return nothing
end

@kernel inbounds = true function soil_carbon_decomposition_kernel!(
    litter::AbstractMatrix{T},
    litter_rate::AbstractVector{T},
    litter_environment::AbstractMatrix{T},
    decomposed_litter::AbstractMatrix{T},
    fast::AbstractMatrix{T},
    slow::AbstractMatrix{T},
    soil_environment::AbstractMatrix{T},
    decomposed_fast::AbstractMatrix{T},
    decomposed_slow::AbstractMatrix{T},
    shift_fast::AbstractMatrix{T},
    shift_slow::AbstractMatrix{T},
    litter_to_fast::AbstractMatrix{T},
    litter_to_slow::AbstractMatrix{T},
    heterotrophic_respiration::AbstractVector{T},
    atmospheric_fraction::T,
    fast_fraction::T,
    fast_rate::T,
    slow_rate::T,
    soil_layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    litter_flux = zero(T)
    for pool in 1:3
        decomposition = (one(T) - exp(
            -litter_rate[pool] * litter_environment[pool, cell],
        )) * litter[pool, cell]
        decomposed_litter[pool, cell] = decomposition
        litter[pool, cell] -= decomposition
        litter_flux += decomposition
    end

    respiration = litter_flux * atmospheric_fraction
    for layer in 1:soil_layers
        fast_decomposition = max(
            zero(T),
            (one(T) - exp(-fast_rate * soil_environment[layer, cell])) *
            fast[layer, cell],
        )
        slow_decomposition = max(
            zero(T),
            (one(T) - exp(-slow_rate * soil_environment[layer, cell])) *
            slow[layer, cell],
        )
        to_fast = shift_fast[layer, cell] * litter_flux * fast_fraction *
            (one(T) - atmospheric_fraction)
        to_slow = shift_slow[layer, cell] * litter_flux *
            (one(T) - fast_fraction) * (one(T) - atmospheric_fraction)
        decomposed_fast[layer, cell] = fast_decomposition
        decomposed_slow[layer, cell] = slow_decomposition
        litter_to_fast[layer, cell] = to_fast
        litter_to_slow[layer, cell] = to_slow
        fast[layer, cell] += to_fast - fast_decomposition
        slow[layer, cell] += to_slow - slow_decomposition
        respiration += fast_decomposition + slow_decomposition
    end
    heterotrophic_respiration[cell] = respiration
end
