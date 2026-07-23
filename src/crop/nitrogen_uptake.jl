# Crop root nitrogen uptake kinetics (LPJmL `nuptake_crop`). Per soil layer, potential uptake of a
# mineral-N pool follows Michaelis-Menten saturation with a baseline term, scaled by a root factor
# that includes a parabolic soil-temperature response. This is the tested scalar kinetics; the
# per-layer, root-weighted coupling to the soil mineral-N pools and the demand-limited two-pass
# accounting are assembled in the crop nitrogen coupling (plan Phase 5).

"""
    $(TYPEDSIGNATURES)

Parabolic soil-temperature response of nitrogen uptake (LPJmL): zero at and below `T_0`, peaking at
the optimum, normalized to 1 at the reference temperature `T_r`. `T_m` is the temperature of maximum
response.
"""
@inline function nitrogen_uptake_temperature_response(T_soil::NF, T_0::NF, T_m::NF, T_r::NF) where {NF}
    numerator = (T_soil - T_0) * (NF(2) * T_m - T_0 - T_soil)
    denominator = (T_r - T_0) * (NF(2) * T_m - T_0 - T_r)
    return max(numerator / denominator, zero(NF))
end

"""
    $(TYPEDEF)

Michaelis-Menten kinetics for root uptake of a mineral-nitrogen pool (NO₃ or NH₄).

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrogenUptakeKinetics{NF}
    "Maximum uptake per unit fine-root carbon"
    @param vmax::NF = 1.5 (bounds = Positive,)
    "Saturation-independent baseline uptake term"
    @param kmin::NF = 0.05 (bounds = Nonnegative,)
    "Half-saturation concentration"
    @param Km::NF = 0.70 (bounds = Positive,)
end

CropNitrogenUptakeKinetics(::Type{NF}; kwargs...) where {NF} = CropNitrogenUptakeKinetics{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Potential root uptake of a mineral-N pool with `available` nitrogen (gN/m²), given the pool's
Michaelis-Menten `saturation_scale` (= wsat·layer_depth, carrying the soil water/geometry) and the
`root_factor` (temperature response × root carbon × root fraction). The result is capped at the
available nitrogen and is zero for an empty pool.
"""
@inline function root_nitrogen_uptake_potential(
        k::CropNitrogenUptakeKinetics{NF}, available::NF, saturation_scale::NF, root_factor::NF,
    ) where {NF}
    pool = max(zero(NF), available)
    saturation = pool / (pool + k.Km * saturation_scale)
    potential = k.vmax * (k.kmin + saturation) * root_factor
    return ifelse(pool > zero(NF), min(potential, pool), zero(NF))
end
