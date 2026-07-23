# Crop phenology: the LPJmL heat-unit-driven leaf-area-index trajectory. The leaf area index is
# `LAI = f(fphu)·laimax`, where `fphu ∈ [0,1]` is the fraction of accumulated phenological heat
# units and `f` is the LPJmL LAI-fraction curve (`src/processes/crop/phenology.jl`): a logistic-like
# rise up to the senescence onset `fphusen`, then a power-law decline to a harvest floor.
#
# Continuous-time note: `fphu` is supplied here as the input `phenology_heat_unit_fraction`. Its
# prognostic accumulation `d(heat units)/dt = max(0, T − T_base)` (with vernalization/photoperiod
# modifiers) is integrated in the crop vegetation model (plan Phase 5); this process contributes the
# faithful LAI *shape* of the growing season.

"""
    $(TYPEDEF)

Crop phenology producing the leaf area index from the phenological-heat-unit fraction via the LPJmL
LAI-fraction trajectory. Defaults are the LPJmL temperate-cereals (CFT 1) values.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropPhenology{NF} <: Terrarium.AbstractPhenology{NF}
    "Heat-unit fraction ending the initial LAI phase"
    @param fphuc::NF = 0.05 (bounds = UnitInterval,)
    "LAI fraction at the end of the initial phase"
    @param flaimaxc::NF = 0.05 (bounds = UnitInterval,)
    "Heat-unit fraction at the LAI sigmoid inflection"
    @param fphuk::NF = 0.45 (bounds = UnitInterval,)
    "LAI fraction at the sigmoid inflection"
    @param flaimaxk::NF = 0.95 (bounds = UnitInterval,)
    "Heat-unit fraction at the onset of senescence"
    @param fphusen::NF = 0.70 (bounds = UnitInterval,)
    "LAI fraction retained at harvest"
    @param flaimaxharvest::NF = 0.0 (bounds = UnitInterval,)
    "Senescence-curve shape parameter"
    @param shapesenescencenorm::NF = 2.0 (bounds = Positive,)
    "Maximum leaf area index"
    @param laimax::NF = 7.0 (bounds = Positive,)
end

CropPhenology(::Type{NF}; kwargs...) where {NF} = CropPhenology{NF}(; kwargs...)

Terrarium.variables(::CropPhenology{NF}) where {NF} = (
    Terrarium.auxiliary(:phenology_factor, XY()),
    Terrarium.auxiliary(:leaf_area_index, XY()),
    Terrarium.input(:phenology_heat_unit_fraction, XY(), default = zero(NF)),
)

# ---- scalar primitives (Level III) --------------------------------------------------------

"""
    $(TYPEDSIGNATURES)

LPJmL leaf-area-index fraction as a function of the heat-unit fraction `fphu ∈ [0,1]`: a
logistic-like rise before `fphusen`, then a power-law senescence decline to `flaimaxharvest`.
"""
@inline function compute_lai_fraction(p::CropPhenology{NF}, fphu::NF) where {NF}
    c = p.fphuc / p.flaimaxc - p.fphuc
    k = p.fphuk / p.flaimaxk - p.fphuk
    growth = fphu / (fphu + c * (c / k)^((p.fphuc - fphu) / (p.fphuk - p.fphuc)))
    senescence = ((one(NF) - fphu) / (one(NF) - p.fphusen))^p.shapesenescencenorm *
        (one(NF) - p.flaimaxharvest) + p.flaimaxharvest
    return ifelse(fphu < p.fphusen, growth, senescence)
end

"""$(TYPEDSIGNATURES) Leaf area index from the heat-unit fraction: `f(fphu)·laimax`."""
@inline function compute_crop_lai(p::CropPhenology{NF}, fphu::NF) where {NF}
    return compute_lai_fraction(p, fphu) * p.laimax
end

# ---- kernel functions (Level II) ----------------------------------------------------------

"""$(TYPEDSIGNATURES) Phenology factor (`fphu`) and leaf area index at grid point `(i, j)`."""
Base.@propagate_inbounds function compute_phenology(i, j, grid, fields, p::CropPhenology{NF}) where {NF}
    fphu = clamp(fields.phenology_heat_unit_fraction[i, j], zero(NF), one(NF))
    LAI = compute_crop_lai(p, fphu)
    return fphu, LAI
end

"""$(TYPEDSIGNATURES) Store [`compute_phenology`](@ref) outputs in `out`."""
Base.@propagate_inbounds function compute_phenology!(out, i, j, grid, fields, p::CropPhenology)
    fphu, LAI = compute_phenology(i, j, grid, fields, p)
    out.phenology_factor[i, j, 1] = fphu
    out.leaf_area_index[i, j, 1] = LAI
    return out
end

# ---- interface methods (Level I) ----------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, p::CropPhenology)
    out = Terrarium.auxiliary_fields(state, p)
    fields = get_fields(state, p; except = out)
    launch!(grid, XY, compute_crop_phenology_kernel!, out, fields, p)
    return nothing
end

@kernel inbounds = true function compute_crop_phenology_kernel!(out, grid, fields, p::CropPhenology)
    i, j = @index(Global, NTuple)
    compute_phenology!(out, i, j, grid, fields, p)
end
