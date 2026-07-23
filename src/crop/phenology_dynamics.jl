# Crop phenological-heat-unit accumulation. The prognostic `phenological_heat_units` (°C·day)
# integrates growing degree-days above a base temperature; its ratio to the crop's heat-unit
# requirement is the fraction `fphu ∈ [0,1]` that drives the crop LAI trajectory (see
# `phenology.jl`). This is the continuous-time replacement for LPJmL's daily heat-unit sum:
# d(HU)/dt = max(0, T_air − T_base), applied per second so integration over a day adds the daily
# growing degree-days. (Vernalization/photoperiod modifiers are future work.)

"""
    $(TYPEDEF)

Prognostic crop phenological-heat-unit accumulation.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropPhenologyDynamics{NF} <: Terrarium.AbstractProcess{NF}
    "Total phenological heat units (°C·day) required to reach maturity"
    @param heat_unit_requirement::NF = 1400.0 (bounds = Positive,)
    "Base temperature above which heat units accumulate"
    @param base_temperature::NF = 0.0 (units = u"°C",)
end

CropPhenologyDynamics(::Type{NF}; kwargs...) where {NF} = CropPhenologyDynamics{NF}(; kwargs...)

Terrarium.variables(::CropPhenologyDynamics{NF}) where {NF} = (
    Terrarium.prognostic(:phenological_heat_units, XY()),
    Terrarium.auxiliary(:phenology_heat_unit_fraction, XY()),
    Terrarium.input(:air_temperature, XY(), default = NF(10), units = u"°C"),
)

"""$(TYPEDSIGNATURES) Heat-unit fraction `fphu = clamp(HU / requirement, 0, 1)`."""
@inline function heat_unit_fraction(pd::CropPhenologyDynamics{NF}, heat_units::NF) where {NF}
    return clamp(heat_units / pd.heat_unit_requirement, zero(NF), one(NF))
end

"""$(TYPEDSIGNATURES) Heat-unit accumulation rate (°C·day per second): `max(0, T_air − T_base)/seconds_per_day`."""
@inline function heat_unit_rate(pd::CropPhenologyDynamics{NF}, air_temperature::NF) where {NF}
    return max(zero(NF), air_temperature - pd.base_temperature) / Terrarium.seconds_per_day(NF)
end

# ---- interface methods --------------------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, pd::CropPhenologyDynamics, args...)
    out = Terrarium.auxiliary_fields(state, pd)
    fields = get_fields(state, pd; except = out)
    launch!(grid, XY, compute_heat_unit_fraction_kernel!, out, fields, pd)
    return nothing
end

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_tendencies!(state, grid, pd::CropPhenologyDynamics, args...)
    tend = Terrarium.tendency_fields(state, pd)
    fields = get_fields(state, pd)
    launch!(grid, XY, compute_heat_unit_tendency_kernel!, tend, fields, pd)
    return nothing
end

@kernel inbounds = true function compute_heat_unit_fraction_kernel!(out, grid, fields, pd::CropPhenologyDynamics)
    i, j = @index(Global, NTuple)
    out.phenology_heat_unit_fraction[i, j, 1] = heat_unit_fraction(pd, fields.phenological_heat_units[i, j])
end

@kernel inbounds = true function compute_heat_unit_tendency_kernel!(tend, grid, fields, pd::CropPhenologyDynamics)
    i, j = @index(Global, NTuple)
    tend.phenological_heat_units[i, j, 1] = heat_unit_rate(pd, fields.air_temperature[i, j])
end
