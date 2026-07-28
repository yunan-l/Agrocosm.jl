"""
soil_carbon!(crop, soil)

Update litter and soil carbon pools and heterotrophic respiration terms.
"""

"""Decompose existing litter and SOM carbon without routing new harvest residues."""
function soil_carbon_decomposition!(soil;
                                    lpjmlparams::LPJmLParams = lpjmlparams,
                                    soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_decomp_response!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    T = eltype(soil_carbon_prognostic(soil).litter)
    launch_custom!(
        soil_carbon_decomposition_kernel!,
        soil_carbon_prognostic(soil).litter,
        size(soil_carbon_prognostic(soil).litter, 2),
        soil_carbon_auxiliary(soil).litter_response,
        soil_decomposition_auxiliary(soil).litter_response,
        soil_carbon_fluxes(soil).decomposed_litter,
        soil_carbon_prognostic(soil).fast,
        soil_carbon_prognostic(soil).slow,
        soil_decomposition_auxiliary(soil).response,
        soil_carbon_fluxes(soil).decomposed_fast,
        soil_carbon_fluxes(soil).decomposed_slow,
        soil_decomposition_input(soil).shift_fast,
        soil_decomposition_input(soil).shift_slow,
        soil_carbon_fluxes(soil).litter_to_fast,
        soil_carbon_fluxes(soil).litter_to_slow,
        soil_carbon_fluxes(soil).heterotrophic_respiration,
        T(lpjmlparams.atmfrac),
        T(lpjmlparams.fastfrac),
        T(lpjmlparams.k_soil10.fast),
        T(lpjmlparams.k_soil10.slow),
        size(soil_carbon_prognostic(soil).fast, 1),
    )
    return nothing
end

"""
    soil_carbon!(crop, soil; ...)

Compatibility entry point for the former combined operation: decompose the
carbon pools, then route harvest-day residues without decomposing those new
residues until the following day.
"""
function soil_carbon!(crop,
                      soil;
                      lpjmlparams::LPJmLParams = lpjmlparams,
                      soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_carbon_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    route_harvest_carbon_input!(soil, crop)
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
        decomposition = -expm1(
            -litter_rate[pool] * litter_environment[pool, cell],
        ) * litter[pool, cell]
        decomposed_litter[pool, cell] = decomposition
        litter[pool, cell] -= decomposition
        litter_flux += decomposition
    end

    respiration = litter_flux * atmospheric_fraction
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
        respiration += fast_decomposition + slow_decomposition
    end
    heterotrophic_respiration[cell] = respiration
end
