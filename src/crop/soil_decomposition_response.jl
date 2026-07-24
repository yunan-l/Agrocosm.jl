# Shared soil-decomposition response (LPJmL `soil_decomp_response`): the temperature and moisture
# rate modifiers that scale heterotrophic decomposition of the soil carbon/nitrogen pools. The
# combined response is a Lloyd-Taylor temperature function (normalized to 1 at 10 °C) times a cubic
# soil-moisture polynomial, clamped to [0,1]. This is the tested scalar physics; it is applied to the
# soil C-N pools in the crop soil-biogeochemistry coupling.

"""
    $(TYPEDEF)

Parameters of the LPJmL soil-decomposition temperature/moisture response.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropSoilDecompositionResponse{NF}
    "Lloyd-Taylor activation parameter (K)"
    @param e0::NF = 308.56 (bounds = Positive,)
    "Lloyd-Taylor reference-temperature offset (K); response is 1 at 10 °C"
    @param temp_response::NF = 56.02 (bounds = Positive,)
    "Intercept of the soil-moisture response polynomial"
    @param intercept::NF = 0.04021601
    "Linear soil-moisture response coefficient"
    @param moist1::NF = 0.71890122
    "Quadratic soil-moisture response coefficient"
    @param moist2::NF = 4.26937932
    "Cubic soil-moisture response coefficient"
    @param moist3::NF = -5.00505434
end

CropSoilDecompositionResponse(::Type{NF}; kwargs...) where {NF} = CropSoilDecompositionResponse{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Lloyd-Taylor temperature response `exp(e0·(1/temp_response − 1/(T + temp_response − 10)))` at air/soil
temperature `T` (°C); equals 1 at the 10 °C reference.
"""
@inline function soil_decomposition_temperature_response(r::CropSoilDecompositionResponse{NF}, temperature::NF) where {NF}
    return exp(r.e0 * (one(NF) / r.temp_response - one(NF) / (temperature + r.temp_response - NF(10))))
end

"""$(TYPEDSIGNATURES) Cubic soil-moisture response polynomial at relative moisture `m ∈ [0,1]`."""
@inline function soil_decomposition_moisture_response(r::CropSoilDecompositionResponse{NF}, moisture::NF) where {NF}
    return r.intercept + r.moist1 * moisture + r.moist2 * moisture^2 + r.moist3 * moisture^3
end

"""$(TYPEDSIGNATURES) Combined decomposition rate modifier (temperature × moisture), clamped to [0,1]."""
@inline function soil_decomposition_response(r::CropSoilDecompositionResponse{NF}, temperature::NF, moisture::NF) where {NF}
    combined = soil_decomposition_temperature_response(r, temperature) * soil_decomposition_moisture_response(r, moisture)
    return clamp(combined, zero(NF), one(NF))
end
