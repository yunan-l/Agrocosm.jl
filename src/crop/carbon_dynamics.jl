# Crop vegetation carbon dynamics: a prognostic vegetation carbon pool and the balanced leaf
# area index derived from it. Follows the PALADYN structure (Willeit & Ganopolski 2016) that
# Terrarium's VegetationCarbon integrates, with two crop-oriented corrections:
#
#   1. Turnover rates are applied in dimensionally-consistent per-second units (Terrarium's
#      PALADYNCarbonDynamics carries them as yr⁻¹ but the timestepper integrates per-second — a
#      flagged inconsistency that drives carbon, and hence LAI, negative).
#   2. The balanced LAI is clamped to be non-negative.
#
# This is the coupled-model-enabling carbon→LAI mapping (leaf carbon × specific leaf area). The
# fully LPJmL-faithful crop LAI additionally caps this by the phenological heat-unit trajectory
# (flaimax·laimax); that requires the crop phenology's prognostic heat-unit accumulation and is
# wired in the crop vegetation model.

"""
    $(TYPEDEF)

Crop vegetation carbon dynamics. Prognostic `carbon_vegetation`; the balanced leaf area index is
`LAI_b = C_veg / (2/SLA + awl)` (clamped ≥ 0), and `d C_veg/dt = (1 - λ_NPP)·NPP - Λ_loc` with
turnover applied per second.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropCarbonDynamics{NF} <: Terrarium.AbstractVegetationCarbonDynamics{NF}
    "Specific leaf area"
    @param SLA::NF = 10.0 (units = u"m^2/kg", bounds = Positive)
    "Allometric coefficient relating stem carbon to leaf area"
    @param awl::NF = 2.0 (units = u"kg/m^2", bounds = Positive)
    "Minimum leaf area index for NPP partitioning"
    @param LAI_min::NF = 1.0 (bounds = Positive,)
    "Maximum leaf area index for NPP partitioning"
    @param LAI_max::NF = 6.0 (bounds = Positive,)
    "Leaf turnover rate"
    @param γL::NF = 0.3 (units = u"yr^-1", bounds = Positive)
    "Root turnover rate"
    @param γR::NF = 0.3 (units = u"yr^-1", bounds = Positive)
    "Stem turnover rate"
    @param γS::NF = 0.05 (units = u"yr^-1", bounds = Positive)
end

CropCarbonDynamics(::Type{NF}; kwargs...) where {NF} = CropCarbonDynamics{NF}(; kwargs...)

Terrarium.variables(::CropCarbonDynamics) = (
    Terrarium.prognostic(:carbon_vegetation, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:balanced_leaf_area_index, XY()),
    Terrarium.input(:net_primary_production, XY(), units = u"kg/m^2/s"),
)

# Seconds per (LPJmL fixed-365) year, for turnover-rate unit conversion.
@inline seconds_per_year(::Type{NF}) where {NF} = NF(365) * Terrarium.seconds_per_day(NF)

# ---- scalar primitives (Level III) --------------------------------------------------------

"""$(TYPEDSIGNATURES) Balanced leaf area index from the vegetation carbon pool (clamped ≥ 0)."""
@inline function compute_balanced_lai(carbon::CropCarbonDynamics{NF}, C_veg::NF) where {NF}
    return max(zero(NF), C_veg / (NF(2) / carbon.SLA + carbon.awl))
end

"""$(TYPEDSIGNATURES) NPP partitioning factor λ_NPP between growth and spreading."""
@inline function compute_lambda_npp(carbon::CropCarbonDynamics{NF}, LAI_b::NF) where {NF}
    return ifelse(
        LAI_b < carbon.LAI_min, zero(NF),
        ifelse(LAI_b ≤ carbon.LAI_max, (LAI_b - carbon.LAI_min) / (carbon.LAI_max - carbon.LAI_min), one(NF)),
    )
end

"""$(TYPEDSIGNATURES) Local litterfall rate Λ_loc (kgC/m²/s) from per-second turnover rates."""
@inline function compute_litterfall(carbon::CropCarbonDynamics{NF}, LAI_b::NF) where {NF}
    per_second = one(NF) / seconds_per_year(NF)
    return (carbon.γL / carbon.SLA + carbon.γR / carbon.SLA + carbon.γS * carbon.awl) * per_second * LAI_b
end

"""$(TYPEDSIGNATURES) Vegetation carbon tendency `d C_veg/dt` (kgC/m²/s)."""
@inline function compute_carbon_tendency(carbon::CropCarbonDynamics{NF}, LAI_b::NF, NPP::NF) where {NF}
    λ_NPP = compute_lambda_npp(carbon, LAI_b)
    Λ_loc = compute_litterfall(carbon, LAI_b)
    return (one(NF) - λ_NPP) * NPP - Λ_loc
end

# ---- kernel functions (Level II) ----------------------------------------------------------

Base.@propagate_inbounds function compute_veg_carbon_auxiliary!(out, i, j, grid, fields, carbon::CropCarbonDynamics)
    out.balanced_leaf_area_index[i, j, 1] = compute_balanced_lai(carbon, fields.carbon_vegetation[i, j])
    return nothing
end

Base.@propagate_inbounds function compute_veg_carbon_tendencies!(tend, i, j, grid, fields, carbon::CropCarbonDynamics)
    LAI_b = fields.balanced_leaf_area_index[i, j]
    NPP = fields.net_primary_production[i, j]
    tend.carbon_vegetation[i, j, 1] = compute_carbon_tendency(carbon, LAI_b, NPP)
    return nothing
end

# ---- interface methods (Level I) ----------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, carbon::CropCarbonDynamics, args...)
    out = Terrarium.auxiliary_fields(state, carbon)
    fields = get_fields(state, carbon; except = out)
    launch!(grid, XY, compute_crop_carbon_auxiliary_kernel!, out, fields, carbon)
    return nothing
end

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_tendencies!(state, grid, carbon::CropCarbonDynamics, args...)
    out = Terrarium.tendency_fields(state, carbon)
    fields = get_fields(state, carbon)
    launch!(grid, XY, compute_crop_carbon_tendencies_kernel!, out, fields, carbon)
    return nothing
end

@kernel inbounds = true function compute_crop_carbon_auxiliary_kernel!(out, grid, fields, carbon::CropCarbonDynamics)
    i, j = @index(Global, NTuple)
    compute_veg_carbon_auxiliary!(out, i, j, grid, fields, carbon)
end

@kernel inbounds = true function compute_crop_carbon_tendencies_kernel!(tend, grid, fields, carbon::CropCarbonDynamics)
    i, j = @index(Global, NTuple)
    compute_veg_carbon_tendencies!(tend, i, j, grid, fields, carbon)
end
