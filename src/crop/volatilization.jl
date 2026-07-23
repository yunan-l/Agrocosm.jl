# Crop soil ammonia volatilization (LPJmL `nitrogen_transform`). The NH₃ flux from the top soil layer
# is the aqueous NH₃ concentration (set by soil pH and the temperature-dependent dissociation) times
# Henry's constant and a wind-driven mass-transfer coefficient, capped at the top-layer ammonium.
# Tested scalar physics; the multi-layer accounting is assembled in the crop soil biogeochemistry
# (plan Phase 5).

"""
    $(TYPEDEF)

Crop soil ammonia-volatilization parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropVolatilization{NF}
    "Characteristic length scale for the mass-transfer coefficient"
    @param length_scale::NF = 1.0 (units = u"m", bounds = Positive)
end

CropVolatilization(::Type{NF}; kwargs...) where {NF} = CropVolatilization{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Ammonia volatilization flux from the top soil layer, given the air temperature (°C), wind speed,
soil pH, top-layer ammonium, and top-layer depth. Capped at the available ammonium.
"""
@inline function ammonia_volatilization(v::CropVolatilization{NF}, air_temperature::NF, wind_speed::NF, soil_ph::NF, ammonium_top::NF, layer_depth::NF) where {NF}
    kelvin = air_temperature + NF(273.15)
    ammonium = max(zero(NF), ammonium_top)
    dissociation = NF(10)^(NF(0.05) - NF(2788) / kelvin)
    aqueous_fraction = one(NF) / (one(NF) + NF(10)^(-soil_ph) / max(dissociation, eps(NF)))
    aqueous_nh3 = aqueous_fraction * ammonium / max(layer_depth, eps(NF)) * NF(1000)
    henry = NF(0.2138) / kelvin * NF(10)^(NF(6.123) - NF(1825) / kelvin)
    mass_transfer = NF(0.000612) * max(zero(NF), wind_speed)^NF(0.8) *
        kelvin^NF(0.382) * v.length_scale^NF(-0.2)
    return clamp(NF(86400) * mass_transfer * henry * aqueous_nh3, zero(NF), ammonium)
end
