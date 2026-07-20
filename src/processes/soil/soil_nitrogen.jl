"""
soil_nitrogen!(crop_cal, crop, soil)

Update litter and soil nitrogen pools and crop-available mineral nitrogen.
"""
function soil_nitrogen_reference!(crop_cal::CropCalendar,
                                  soil::Soil;
                                  air_temperature = nothing,
                                  wind_speed = nothing,
                                  lpjmlparams::LPJmLParams = lpjmlparams,
                                  soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack atmfrac, fastfrac, k_soil10 = lpjmlparams
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # compute soil carbon: litter carbon and soil carbon
    # soil.nitrogen.decomposed_litter = (1.0f0 .- exp.(-soil.nitrogen.litter_response / 100)) .* soil.nitrogen.litter

    soil.nitrogen.decomposed_litter .=
        -expm1.(-soil.nitrogen.litter_response .* soil.decomposition.litter_response) .* soil.nitrogen.litter
    soil.nitrogen.litter .-= soil.nitrogen.decomposed_litter

    route_harvest_nitrogen_input_reference!(soil, crop_cal)

    decomposed_litter = soil.decomposition.surface_scratch_1
    @views decomposed_litter .=
        soil.nitrogen.decomposed_litter[1, :] .+
        soil.nitrogen.decomposed_litter[2, :] .+
        soil.nitrogen.decomposed_litter[3, :]
    soil.nitrogen.litter_to_fast .= soil.nitrogen.shift_fast .*
        reshape(decomposed_litter, 1, :) .* fastfrac .* (1.0f0 - atmfrac)
    soil.nitrogen.litter_to_slow .= soil.nitrogen.shift_slow .*
        reshape(decomposed_litter, 1, :) .* (1.0f0 - fastfrac) .* (1.0f0 - atmfrac)

    # soil.nitrogen.decomposed_fast = (1.0f0 .- exp.(-soil.response_fastn .* response / 50)) .* soil.nitrogen.fast
    soil.nitrogen.decomposed_fast .= max.(
        0.0f0,
        -expm1.(-k_soil10.fast .* soil.decomposition.response) .* soil.nitrogen.fast,
    )
    soil.nitrogen.fast .+= soil.nitrogen.litter_to_fast .- soil.nitrogen.decomposed_fast

    # soil.nitrogen.decomposed_slow = (1.0f0 .- exp.(-soil.response_slown .* response / 10)) .* soil.nitrogen.slow
    soil.nitrogen.decomposed_slow .= max.(
        0.0f0,
        -expm1.(-k_soil10.slow .* soil.decomposition.response) .* soil.nitrogen.slow,
    )
    soil.nitrogen.slow .+= soil.nitrogen.litter_to_slow .- soil.nitrogen.decomposed_slow

    nitrogen_transform!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )

end

function soil_nitrogen!(crop_cal::CropCalendar,
                        soil::Soil;
                        air_temperature = nothing,
                        wind_speed = nothing,
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        soil_decomp_params::SoilDecompParams = soil_decomp_params)
    T = eltype(soil.nitrogen.litter)
    launch_custom!(
        soil_nitrogen_decomposition_kernel!,
        soil.nitrogen.litter,
        size(soil.nitrogen.litter, 2),
        soil.nitrogen.litter_response,
        soil.decomposition.litter_response,
        soil.nitrogen.decomposed_litter,
        soil.nitrogen.fast,
        soil.nitrogen.slow,
        soil.decomposition.response,
        soil.nitrogen.decomposed_fast,
        soil.nitrogen.decomposed_slow,
        soil.nitrogen.shift_fast,
        soil.nitrogen.shift_slow,
        soil.nitrogen.litter_to_fast,
        soil.nitrogen.litter_to_slow,
        T(lpjmlparams.atmfrac),
        T(lpjmlparams.fastfrac),
        T(lpjmlparams.k_soil10.fast),
        T(lpjmlparams.k_soil10.slow),
        size(soil.nitrogen.fast, 1),
    )
    route_harvest_nitrogen_input!(soil, crop_cal)
    nitrogen_transform!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function soil_nitrogen_decomposition_kernel!(
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
    atmospheric_fraction::T,
    fast_fraction::T,
    fast_rate::T,
    slow_rate::T,
    soil_layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    litter_flux = zero(T)
    for pool in 1:3
        decomposition = -expm1(
            -litter_rate[pool] * litter_environment[pool, cell],
        ) * litter[pool, cell]
        decomposed_litter[pool, cell] = decomposition
        litter[pool, cell] -= decomposition
        litter_flux += decomposition
    end

    for layer in 1:soil_layers
        fast_decomposition = max(
            zero(T),
            -expm1(-fast_rate * soil_environment[layer, cell]) *
            fast[layer, cell],
        )
        slow_decomposition = max(
            zero(T),
            -expm1(-slow_rate * soil_environment[layer, cell]) *
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
    end
end
