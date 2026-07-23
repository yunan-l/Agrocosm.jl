# Crop soil nitrogen mineralization / immobilization (LPJmL `nitrogen_transform`). Decomposing
# organic matter with a wide C:N ratio immobilizes mineral nitrogen to meet the microbial demand
# `litter_C/soil_CN − litter_N`; the actual immobilization is Michaelis-Menten-limited by the
# available mineral nitrogen concentration and capped at the available pool. (Gross mineralization
# itself is the nitrogen released by the decomposed carbon pools — an input from decomposition.)
# Tested scalar physics; the multi-layer fast/slow accounting is assembled in the crop soil
# biogeochemistry (plan Phase 5).

"""
    $(TYPEDEF)

Crop soil mineralization/immobilization parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrogenMineralization{NF}
    "Soil organic-matter C:N ratio setting the immobilization demand"
    @param soil_cn_ratio::NF = 15.0 (bounds = Positive,)
    "Half-saturation concentration for immobilization"
    @param immobilization_k::NF = 5.0e-3 (bounds = Positive,)
end

CropNitrogenMineralization(::Type{NF}; kwargs...) where {NF} = CropNitrogenMineralization{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Microbial immobilization demand `max(0, litter_carbon/soil_CN − litter_nitrogen)`. Positive when the
decomposing litter is nitrogen-poor relative to the soil C:N ratio; zero when it is nitrogen-rich
(net mineralization).
"""
@inline function immobilization_demand(m::CropNitrogenMineralization{NF}, litter_carbon::NF, litter_nitrogen::NF) where {NF}
    return max(zero(NF), litter_carbon / m.soil_cn_ratio - litter_nitrogen)
end

"""
    $(TYPEDSIGNATURES)

Michaelis-Menten limitation factor on immobilization from the available mineral-nitrogen
concentration `available/layer_depth·1000`: `concentration / (immobilization_k + concentration)`.
"""
@inline function immobilization_limitation(m::CropNitrogenMineralization{NF}, available_nitrogen::NF, layer_depth::NF) where {NF}
    concentration = max(zero(NF), available_nitrogen) / max(layer_depth, eps(NF)) * NF(1000)
    return concentration / (m.immobilization_k + concentration)
end

"""
    $(TYPEDSIGNATURES)

Actual immobilized nitrogen: the demand scaled by the Michaelis-Menten limitation, capped at the
available mineral nitrogen.
"""
@inline function immobilized_nitrogen(m::CropNitrogenMineralization{NF}, demand::NF, available_nitrogen::NF, layer_depth::NF) where {NF}
    limitation = immobilization_limitation(m, available_nitrogen, layer_depth)
    return min(max(zero(NF), demand) * limitation, max(zero(NF), available_nitrogen))
end
