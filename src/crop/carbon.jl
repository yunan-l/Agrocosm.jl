# Prognostic crop carbon pool. Total plant biomass (kgC/m²) accumulates net primary production and is
# partitioned into organs each step using the ported allocation and respiration primitives:
#
#   root_carbon    = root_allocation_fraction(fphu, df)·biomass          (df = 100·β soil-water stress)
#   leaf_carbon    = min(LAI/SLA, biomass − root_carbon)                 (phenology LAI, carbon-limited)
#   storage_carbon = max(0, biomass − root_carbon − leaf_carbon)
#   NPP = GPP − Rm − Rg,   Rm = maintenance respiration,   Rg = growth respiration
#   d(biomass)/dt = NPP
#
# All rates are in kgC/m²/s (GPP is already net of leaf respiration in Terrarium's convention, i.e.
# net leaf-level assimilation × 1e-3); the maintenance respiration is converted from its per-day form.
# This closes the crop carbon loop; the LAI-feedback deficit and multi-pool nitrogen are future work.

"""
    $(TYPEDEF)

Prognostic crop carbon pool and organ partitioning.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropCarbon{NF} <: Terrarium.AbstractProcess{NF}
    "Root carbon-allocation parameters"
    allocation::CropCarbonAllocation{NF} = CropCarbonAllocation(NF)
    "Maintenance-respiration parameters"
    maintenance::CropMaintenanceRespiration{NF} = CropMaintenanceRespiration(NF)
    "Growth-respiration parameters"
    growth::CropGrowthRespiration{NF} = CropGrowthRespiration(NF)
    "Specific leaf area (m² leaf per kgC)"
    specific_leaf_area::NF = 30.0
end

CropCarbon(::Type{NF}; kwargs...) where {NF} = CropCarbon{NF}(; kwargs...)

Terrarium.variables(::CropCarbon{NF}) where {NF} = (
    Terrarium.prognostic(:crop_biomass, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:leaf_carbon, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:root_carbon, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:storage_carbon, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:net_primary_production, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:gross_primary_production, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:leaf_area_index, XY()),
    Terrarium.input(:phenology_heat_unit_fraction, XY()),
    Terrarium.input(:soil_moisture_limiting_factor, XY(), default = NF(1)),
    Terrarium.input(:air_temperature, XY(), default = NF(10), units = u"°C"),
)

"""
    $(TYPEDSIGNATURES)

Partition biomass into `(leaf_carbon, root_carbon, storage_carbon)` and compute net primary
production, given the heat-unit fraction, soil-water stress β, leaf area index, temperature, and GPP.
"""
@inline function crop_carbon_budget(
        c::CropCarbon{NF}, biomass::NF, fphu::NF, β::NF, leaf_area_index::NF, temperature::NF, gpp::NF,
    ) where {NF}
    # Water/nitrogen stress index for allocation on the LPJmL 0–100 scale (low stress → fewer roots).
    df = NF(100) * clamp(β, zero(NF), one(NF))
    root_fraction = root_allocation_fraction(c.allocation, fphu, df)
    root_carbon = max(zero(NF), root_fraction * biomass)
    leaf_carbon = leaf_carbon_from_lai(leaf_area_index, c.specific_leaf_area, biomass - root_carbon)
    storage_carbon = max(zero(NF), biomass - root_carbon - leaf_carbon)
    # Maintenance respiration (per-day form → per second); mobile pool folded into storage here.
    Rm = maintenance_respiration(c.maintenance, root_carbon, storage_carbon, zero(NF), temperature, temperature) /
        Terrarium.seconds_per_day(NF)
    npp = net_primary_production(c.growth, gpp, Rm)
    return leaf_carbon, root_carbon, storage_carbon, npp
end

# ---- interface methods --------------------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, c::CropCarbon, args...)
    out = Terrarium.auxiliary_fields(state, c)
    fields = get_fields(state, c; except = out)
    launch!(grid, XY, compute_crop_carbon_auxiliary_kernel!, out, fields, c)
    return nothing
end

""" $(TYPEDSIGNATURES) Integrate the crop biomass: d(biomass)/dt = NPP. """
function Terrarium.compute_tendencies!(state, grid, c::CropCarbon, args...)
    tend = Terrarium.tendency_fields(state, c)
    fields = get_fields(state, c)
    launch!(grid, XY, compute_crop_carbon_tendency_kernel!, tend, fields, c)
    return nothing
end

@kernel inbounds = true function compute_crop_carbon_auxiliary_kernel!(out, grid, fields, c::CropCarbon)
    i, j = @index(Global, NTuple)
    leaf, root, storage, npp = crop_carbon_budget(
        c, fields.crop_biomass[i, j], fields.phenology_heat_unit_fraction[i, j],
        fields.soil_moisture_limiting_factor[i, j], fields.leaf_area_index[i, j],
        fields.air_temperature[i, j], fields.gross_primary_production[i, j],
    )
    out.leaf_carbon[i, j, 1] = leaf
    out.root_carbon[i, j, 1] = root
    out.storage_carbon[i, j, 1] = storage
    out.net_primary_production[i, j, 1] = npp
end

@kernel inbounds = true function compute_crop_carbon_tendency_kernel!(tend, grid, fields, c::CropCarbon)
    i, j = @index(Global, NTuple)
    tend.crop_biomass[i, j, 1] = fields.net_primary_production[i, j]
end
