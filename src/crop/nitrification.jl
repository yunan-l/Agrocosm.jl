# Crop soil nitrification (LPJmL `nitrogen_transform`). Gross nitrification of ammonium to nitrate is
# the product of a maximum rate, the ammonium stock, and three environmental factors — a peaked
# water-filled-pore-space (WFPS) moisture response, a Gaussian soil-temperature response, and an
# atan soil-pH response — capped at the available ammonium. A fixed fraction is lost as N₂O. Tested
# scalar physics; the multi-layer soil mineral-N accounting is assembled in the crop soil
# biogeochemistry (plan Phase 5).

"""
    $(TYPEDEF)

Crop soil nitrification parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrification{NF}
    "Maximum fraction of soil ammonium nitrified"
    @param k_max::NF = 0.10 (bounds = UnitInterval,)
    "Fraction of gross nitrification lost as N₂O"
    @param k_2::NF = 0.01 (bounds = UnitInterval,)
    "WFPS moisture-response parameter a"
    @param a::NF = 0.45
    "WFPS moisture-response parameter b"
    @param b::NF = 1.27
    "WFPS moisture-response parameter c"
    @param c::NF = 0.0012
    "WFPS moisture-response parameter d"
    @param d::NF = 2.84
end

CropNitrification(::Type{NF}; kwargs...) where {NF} = CropNitrification{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Peaked water-filled-pore-space (WFPS) moisture response for nitrification, zero outside its support.
"""
@inline function nitrification_moisture_factor(n::CropNitrification{NF}, water_filled_pore_space::NF) where {NF}
    n_nit = n.a - n.b
    m_nit = n.a - n.c
    z_nit = n.d * (n.b - n.a) / (n.a - n.c)
    # Clamp each base to ≥ 0 before the (non-integer) power: outside the support one base is
    # negative and `base^d` would throw a DomainError. A zeroed base makes the product 0, exactly
    # the LPJmL "0 outside support" behaviour, while staying throw-free (kernel/Reactant-safe).
    base_1 = max(zero(NF), (water_filled_pore_space - n.b) / n_nit)
    base_2 = max(zero(NF), (water_filled_pore_space - n.c) / m_nit)
    return base_1^z_nit * base_2^n.d
end

"""$(TYPEDSIGNATURES) Gaussian soil-temperature response for nitrification (peak at 18.79 °C)."""
@inline nitrification_temperature_factor(::CropNitrification{NF}, soil_temperature::NF) where {NF} =
    exp(-(soil_temperature - NF(18.79))^2 / NF(2 * 8.26 * 8.26))

"""$(TYPEDSIGNATURES) Soil-pH response for nitrification."""
@inline nitrification_ph_factor(::CropNitrification{NF}, soil_ph::NF) where {NF} =
    NF(0.56) + atan(NF(pi) * NF(0.45) * (soil_ph - NF(5))) / NF(pi)

"""
    $(TYPEDSIGNATURES)

Gross nitrification (NH₄→NO₃) and the associated N₂O loss, given the ammonium stock, WFPS, soil
temperature (°C), and soil pH. Returns `(gross_nitrification, n2o_loss)`, with gross nitrification
capped at the available ammonium.
"""
@inline function gross_nitrification(n::CropNitrification{NF}, ammonium::NF, water_filled_pore_space::NF, soil_temperature::NF, soil_ph::NF) where {NF}
    moisture = nitrification_moisture_factor(n, water_filled_pore_space)
    temperature = nitrification_temperature_factor(n, soil_temperature)
    ph = nitrification_ph_factor(n, soil_ph)
    gross = clamp(n.k_max * ammonium * temperature * moisture * ph, zero(NF), ammonium)
    return gross, n.k_2 * gross
end
