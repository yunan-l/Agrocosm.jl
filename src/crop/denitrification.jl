# Crop soil denitrification (LPJmL `nitrogen_transform`). Gross denitrification of nitrate is the
# product of a water-filled-pore-space moisture factor, a soil-carbon availability factor (driven by
# the fast + slow organic carbon and a soil-temperature response), and the nitrate stock; the gaseous
# loss is split into N₂O and N₂. Tested scalar physics; the multi-layer soil mineral-N accounting is
# assembled in the crop soil biogeochemistry.

"""
    $(TYPEDEF)

Crop soil denitrification parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropDenitrification{NF}
    "Denitrification shape factor on the organic-carbon availability"
    @param CDN::NF = 1.2 (bounds = Positive,)
    "Fraction of denitrified nitrogen emitted as N₂O (rest as N₂)"
    @param n2o_fraction::NF = 0.11 (bounds = UnitInterval,)
end

CropDenitrification(::Type{NF}; kwargs...) where {NF} = CropDenitrification{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Soil-temperature response for denitrification: zero above 45.9 °C, a peaked polynomial for
`0 < T ≤ 45.9`, and a constant `0.0326` at or below 0 °C. The powers use a non-negative base so the
function is throw-free (kernel/Reactant-safe).
"""
@inline function denitrification_temperature_factor(::CropDenitrification{NF}, soil_temperature::NF) where {NF}
    t = max(zero(NF), soil_temperature)
    active = max(zero(NF), NF(0.0326) + NF(0.00351) * t^NF(1.652) - (t / NF(41.748))^NF(7.19))
    return ifelse(soil_temperature > NF(45.9), zero(NF), ifelse(soil_temperature > zero(NF), active, NF(0.0326)))
end

"""$(TYPEDSIGNATURES) Water-filled-pore-space moisture factor for denitrification (capped at 1)."""
@inline denitrification_moisture_factor(::CropDenitrification{NF}, water_filled_pore_space::NF) where {NF} =
    min(one(NF), NF(6.664096e-10) * exp(NF(20.92912) * water_filled_pore_space))

"""
    $(TYPEDSIGNATURES)

Gross denitrification and its N₂O/N₂ split, given the nitrate stock, soil temperature (°C),
water-filled pore space, and the organic (fast + slow) carbon. Returns
`(gross_denitrification, n2o_loss, n2_loss)`, capped at the available nitrate.
"""
@inline function gross_denitrification(dn::CropDenitrification{NF}, nitrate::NF, soil_temperature::NF, water_filled_pore_space::NF, organic_carbon::NF) where {NF}
    temperature = denitrification_temperature_factor(dn, soil_temperature)
    moisture = denitrification_moisture_factor(dn, water_filled_pore_space)
    carbon_factor = max(zero(NF), one(NF) - exp(-dn.CDN * temperature * max(zero(NF), organic_carbon)))
    gross = clamp(moisture * carbon_factor * nitrate, zero(NF), nitrate)
    n2o = dn.n2o_fraction * gross
    return gross, n2o, gross - n2o
end
