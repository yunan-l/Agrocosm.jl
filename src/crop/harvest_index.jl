# Crop harvest index (LPJmL `carbon_allocation`). The fraction of accumulated carbon partitioned to
# the harvested (storage) organ: a phenology-driven optimum `fhiopt` (a sigmoid in the heat-unit
# fraction `fphu`) scaled between the minimum (`himin`) and optimal (`hiopt`) harvest index, then
# reduced by the water-deficit sufficiency factor `wdf`. Values > 1 encode above-ground vs total
# biomass ratios for some crops (e.g. sugarcane). Tested scalar physics; applied in the crop carbon
# allocation coupling (plan Phase 5).

"""
    $(TYPEDEF)

Crop harvest-index parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropHarvestIndex{NF}
    "Optimal harvest index (well-watered)"
    @param hiopt::NF = 0.5 (bounds = Positive,)
    "Minimum harvest index under water stress"
    @param himin::NF = 0.2 (bounds = Positive,)
end

CropHarvestIndex(::Type{NF}; kwargs...) where {NF} = CropHarvestIndex{NF}(; kwargs...)

# HI values above 1 (above-ground vs total ratios) are scaled about 1 rather than about 0.
@inline _scale_hi(fhiopt::NF, hi::NF) where {NF} = ifelse(hi > one(NF), fhiopt * (hi - one(NF)) + one(NF), fhiopt * hi)

"""
    $(TYPEDSIGNATURES)

Harvest index from the heat-unit fraction `fphu ∈ [0,1]` and the water-deficit sufficiency factor
`wdf` (high = well-watered): approaches the phenology-scaled optimum at large `wdf` and the scaled
minimum at `wdf = 0`.
"""
@inline function crop_harvest_index(h::CropHarvestIndex{NF}, fphu::NF, wdf::NF) where {NF}
    fhiopt = NF(100) * fphu / (NF(100) * fphu + exp(NF(11.1) - NF(10) * fphu))
    hi_opt = _scale_hi(fhiopt, h.hiopt)
    hi_min = _scale_hi(fhiopt, h.himin)
    water_factor = wdf / (wdf + exp(NF(6.13) - NF(0.0883) * wdf))
    return (hi_opt - hi_min) * water_factor + hi_min
end
