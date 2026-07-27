# Crop leaf-nitrogen limitation of the Rubisco capacity (LPJmL `limit_vcmax_by_nitrogen`).
# Structural leaf nitrogen (ncleaf_min·leaf_carbon) is protected; only the excess supports Rubisco
# activity, which — temperature-scaled — caps the potential Vcmax. This is applied to the crop
# photosynthesis Vcmax within the crop nitrogen coupling; here it is the tested
# scalar physics.

"""
    $(TYPEDEF)

Parameters of the crop leaf-nitrogen limitation of Rubisco capacity.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrogenVcmaxLimit{NF}
    "Minimum (structural) leaf nitrogen-to-carbon ratio protected from Rubisco use"
    @param ncleaf_min::NF = 1 / 58.8 (bounds = Positive,)
    "Temperature sensitivity of the nitrogen requirement for Rubisco activity"
    @param k_temp::NF = 0.0693 (bounds = Positive,)
    "Atmospheric-pressure/psychrometric scaling constant (LPJmL `p`)"
    @param pressure_scale::NF = 25.0 (bounds = Positive,)
end

CropNitrogenVcmaxLimit(::Type{NF}; kwargs...) where {NF} = CropNitrogenVcmaxLimit{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Nitrogen-limited Rubisco capacity and the retained fraction, given the potential Vcmax, the
available leaf nitrogen, leaf carbon, and air temperature (°C). Returns
`(limited_vcmax, nitrogen_limitation ∈ [0,1])`. When there is no potential capacity, nothing is
limited and the retained fraction is zero.
"""
@inline function nitrogen_limited_vcmax(
        limit::CropNitrogenVcmaxLimit{NF}, potential_vcmax::NF,
        available_leaf_nitrogen::NF, leaf_carbon::NF, temperature::NF,
    ) where {NF}
    potential = max(zero(NF), potential_vcmax)
    # Only leaf nitrogen in excess of the structural minimum supports Rubisco.
    rubisco_nitrogen = max(zero(NF), available_leaf_nitrogen - limit.ncleaf_min * leaf_carbon)
    nitrogen_capacity = rubisco_nitrogen /
        exp(-limit.k_temp * (temperature - NF(25))) /
        (limit.pressure_scale * NF(1.0e-3)) * (NF(86400) * NF(12) * NF(1.0e-6))
    capped = min(potential, max(eps(NF), nitrogen_capacity))
    limited = ifelse(potential > zero(NF), capped, potential)
    fraction = ifelse(potential > zero(NF), clamp(limited / max(potential, eps(NF)), zero(NF), one(NF)), zero(NF))
    return limited, fraction
end

"""
    $(TYPEDSIGNATURES)

Nitrogen-supported maximum carboxylation rate (gC/m²/s, matching the photosynthesis `Vc_max`) from the
leaf nitrogen and carbon pools (kgN/m², kgC/m²) and air temperature (°C). This is the Rubisco capacity
the available leaf nitrogen can sustain; capping the light-derived potential `Vc_max` at it applies the
LPJmL leaf-nitrogen Rubisco limitation. Returns a large value when there is no leaf carbon yet (early
growth is not nitrogen-limited). The `1000` converts the revised kg pools to the LPJmL g pools, and the
`1/86400` converts the LPJmL per-day rate to the per-second `Vc_max`.
"""
@inline function nitrogen_supported_vcmax(
        limit::CropNitrogenVcmaxLimit{NF}, leaf_nitrogen::NF, leaf_carbon::NF, temperature::NF,
    ) where {NF}
    leaf_carbon > zero(NF) || return NF(Inf)
    rubisco_nitrogen = max(zero(NF), (leaf_nitrogen - limit.ncleaf_min * leaf_carbon) * NF(1000))
    capacity_per_day = rubisco_nitrogen /
        exp(-limit.k_temp * (temperature - NF(25))) /
        (limit.pressure_scale * NF(1.0e-3)) * (NF(86400) * NF(12) * NF(1.0e-6))
    return capacity_per_day / NF(86400)
end
