# Crop growth respiration and net primary production. A fixed fraction `r_growth` of the assimilate
# remaining after maintenance respiration is spent on growth respiration (LPJmL `r_growth = 0.25`):
# Rg = r_growth·max(0, GPP − Rm), and NPP = GPP − Rm − Rg. This is the tested scalar physics; it is
# applied within the crop carbon coupling.

"""
    $(TYPEDEF)

Crop growth-respiration parameter.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropGrowthRespiration{NF}
    "Fraction of the post-maintenance assimilate consumed by growth respiration"
    @param r_growth::NF = 0.25 (bounds = UnitInterval,)
end

CropGrowthRespiration(::Type{NF}; kwargs...) where {NF} = CropGrowthRespiration{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Growth respiration `Rg = r_growth·max(0, gross_assimilation − maintenance_respiration)`. Only positive
net assimilate incurs growth respiration.
"""
@inline function growth_respiration(g::CropGrowthRespiration{NF}, gross_assimilation::NF, maintenance_respiration::NF) where {NF}
    return g.r_growth * max(zero(NF), gross_assimilation - maintenance_respiration)
end

"""
    $(TYPEDSIGNATURES)

Net primary production `NPP = gross_assimilation − maintenance_respiration − Rg`. Equals
`(1 − r_growth)·(gross − maintenance)` when that is positive, and `gross − maintenance` (a loss) when
negative.
"""
@inline function net_primary_production(g::CropGrowthRespiration{NF}, gross_assimilation::NF, maintenance_respiration::NF) where {NF}
    return gross_assimilation - maintenance_respiration - growth_respiration(g, gross_assimilation, maintenance_respiration)
end
