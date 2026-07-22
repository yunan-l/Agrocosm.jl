"""
soil_nitrogen!(crop, soil)

Update litter and soil nitrogen pools and crop-available mineral nitrogen.
"""
function soil_nitrogen_reference!(crop,
                                  soil;
                                  air_temperature = nothing,
                                  wind_speed = nothing,
                                  lpjmlparams::LPJmLParams = lpjmlparams,
                                  soil_decomp_params::SoilDecompParams = soil_decomp_params
)

    @unpack atmfrac, fastfrac, k_soil10 = lpjmlparams
    @unpack intercept, moist3, moist2, moist1, eps = soil_decomp_params

    # compute soil carbon: litter carbon and soil carbon
    # soil_nitrogen_fluxes(soil).decomposed_litter = (1.0f0 .- exp.(-soil_nitrogen_auxiliary(soil).litter_response / 100)) .* soil_nitrogen_prognostic(soil).litter

    soil_nitrogen_fluxes(soil).decomposed_litter .=
        -expm1.(-soil_nitrogen_auxiliary(soil).litter_response .* soil_decomposition_auxiliary(soil).litter_response) .* soil_nitrogen_prognostic(soil).litter
    soil_nitrogen_prognostic(soil).litter .-= soil_nitrogen_fluxes(soil).decomposed_litter

    route_harvest_nitrogen_input_reference!(soil, crop)

    decomposed_litter = soil_decomposition_workspace(soil).surface_scratch_1
    @views decomposed_litter .=
        soil_nitrogen_fluxes(soil).decomposed_litter[1, :] .+
        soil_nitrogen_fluxes(soil).decomposed_litter[2, :] .+
        soil_nitrogen_fluxes(soil).decomposed_litter[3, :]
    soil_nitrogen_fluxes(soil).litter_to_fast .= soil_decomposition_input(soil).shift_fast .*
        reshape(decomposed_litter, 1, :) .* fastfrac .* (1.0f0 - atmfrac)
    soil_nitrogen_fluxes(soil).litter_to_slow .= soil_decomposition_input(soil).shift_slow .*
        reshape(decomposed_litter, 1, :) .* (1.0f0 - fastfrac) .* (1.0f0 - atmfrac)

    # soil_nitrogen_fluxes(soil).decomposed_fast = (1.0f0 .- exp.(-soil.response_fastn .* response / 50)) .* soil_nitrogen_prognostic(soil).fast
    soil_nitrogen_fluxes(soil).decomposed_fast .= max.(
        0.0f0,
        -expm1.(-k_soil10.fast .* soil_decomposition_auxiliary(soil).response) .* soil_nitrogen_prognostic(soil).fast,
    )
    soil_nitrogen_prognostic(soil).fast .+= soil_nitrogen_fluxes(soil).litter_to_fast .- soil_nitrogen_fluxes(soil).decomposed_fast

    # soil_nitrogen_fluxes(soil).decomposed_slow = (1.0f0 .- exp.(-soil.response_slown .* response / 10)) .* soil_nitrogen_prognostic(soil).slow
    soil_nitrogen_fluxes(soil).decomposed_slow .= max.(
        0.0f0,
        -expm1.(-k_soil10.slow .* soil_decomposition_auxiliary(soil).response) .* soil_nitrogen_prognostic(soil).slow,
    )
    soil_nitrogen_prognostic(soil).slow .+= soil_nitrogen_fluxes(soil).litter_to_slow .- soil_nitrogen_fluxes(soil).decomposed_slow

    nitrogen_transform!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )

end

"""Decompose existing litter and SOM nitrogen without mineral transformations."""
function soil_nitrogen_decomposition!(soil;
                                      lpjmlparams::LPJmLParams = lpjmlparams,
                                      soil_decomp_params::SoilDecompParams = soil_decomp_params,
                                      litter_rate = soil_nitrogen_auxiliary(soil).litter_response,
                                      shift_fast = soil_decomposition_input(soil).shift_fast,
                                      shift_slow = soil_decomposition_input(soil).shift_slow)
    T = eltype(soil_nitrogen_prognostic(soil).litter)
    launch_custom!(
        soil_nitrogen_decomposition_kernel!,
        soil_nitrogen_prognostic(soil).litter,
        size(soil_nitrogen_prognostic(soil).litter, 2),
        litter_rate,
        soil_decomposition_auxiliary(soil).litter_response,
        soil_nitrogen_fluxes(soil).decomposed_litter,
        soil_nitrogen_prognostic(soil).fast,
        soil_nitrogen_prognostic(soil).slow,
        soil_decomposition_auxiliary(soil).response,
        soil_nitrogen_fluxes(soil).decomposed_fast,
        soil_nitrogen_fluxes(soil).decomposed_slow,
        shift_fast,
        shift_slow,
        soil_nitrogen_fluxes(soil).litter_to_fast,
        soil_nitrogen_fluxes(soil).litter_to_slow,
        T(lpjmlparams.atmfrac),
        T(lpjmlparams.fastfrac),
        T(lpjmlparams.k_soil10.fast),
        T(lpjmlparams.k_soil10.slow),
        size(soil_nitrogen_prognostic(soil).fast, 1),
    )
    return nothing
end

"""
    soil_cn_decomposition!(soil; ...)

Execute LPJmL's coupled pre-crop soil stage: compute one shared environmental
response, decompose C and N with identical pool-specific decay fractions,
then mineralize/immobilize and nitrify mineral nitrogen.
"""
function soil_cn_decomposition!(soil;
                                lpjmlparams::LPJmLParams = lpjmlparams,
                                soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_carbon_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    soil_nitrogen_decomposition!(
        soil;
        lpjmlparams = lpjmlparams,
        soil_decomp_params = soil_decomp_params,
        litter_rate = soil_carbon_auxiliary(soil).litter_response,
        shift_fast = soil_decomposition_input(soil).shift_fast,
        shift_slow = soil_decomposition_input(soil).shift_slow,
    )
    mineralize_nitrify!(
        soil;
        lpjmlparams = lpjmlparams,
        shift_fast = soil_decomposition_input(soil).shift_fast,
        shift_slow = soil_decomposition_input(soil).shift_slow,
    )
    return nothing
end

"""Route new harvest-day carbon and nitrogen residues together."""
function route_harvest_residues!(soil, crop)
    route_harvest_carbon_input!(soil, crop)
    route_harvest_nitrogen_input!(soil, crop)
    return nothing
end

"""
    soil_nitrogen!(crop, soil; ...)

Compatibility entry point for the former combined operation.
"""
function soil_nitrogen!(crop,
                        soil;
                        air_temperature = nothing,
                        wind_speed = nothing,
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_nitrogen_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    route_harvest_nitrogen_input!(soil, crop)
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
