# Crop carbon allocation to organs (LPJmL/SWAT `carbon_allocation`). The root fraction of biomass
# starts at `FROOTMAX` and declines through the season (with the heat-unit fraction `fphu`) toward
# `FROOTMAX − FROOTMIN`, modulated by a water/nitrogen stress factor `df` (low stress → more roots).
# Leaf carbon is constrained by the phenological LAI and the specific leaf area. Tested scalar
# physics; the stateful biomass/senescence accounting is assembled in the crop carbon coupling
#.

"""
    $(TYPEDEF)

Crop root carbon-allocation parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropCarbonAllocation{NF}
    "Maximum root fraction of biomass (early season / full stress)"
    @param root_fraction_max::NF = 0.4 (bounds = UnitInterval,)
    "Seasonal reduction of the root fraction"
    @param root_fraction_min::NF = 0.3 (bounds = UnitInterval,)
end

CropCarbonAllocation(::Type{NF}; kwargs...) where {NF} = CropCarbonAllocation{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Root fraction of biomass from the heat-unit fraction `fphu ∈ [0,1]` and the water/nitrogen stress
factor `df` (SWAT-style): `root_fraction_max − root_fraction_min·fphu·df/(df + exp(6.13 − 0.0883·df))`.
At low stress (small `df`) the fraction stays near its maximum (root investment); well-watered late in
the season it approaches `root_fraction_max − root_fraction_min`.
"""
@inline function root_allocation_fraction(a::CropCarbonAllocation{NF}, fphu::NF, df::NF) where {NF}
    stress = df / (df + exp(NF(6.13) - NF(0.0883) * df))
    return a.root_fraction_max - a.root_fraction_min * fphu * stress
end

"""
    $(TYPEDSIGNATURES)

Leaf carbon constrained by the phenological leaf area index and the specific leaf area: the smaller of
`leaf_area_index / specific_leaf_area` and the carbon available for leaves.
"""
@inline function leaf_carbon_from_lai(leaf_area_index::NF, specific_leaf_area::NF, available_carbon::NF) where {NF}
    return min(leaf_area_index / specific_leaf_area, max(zero(NF), available_carbon))
end
