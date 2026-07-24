# Prognostic crop nitrogen pool. Total plant nitrogen (kgN/m²) is acquired in proportion to the net
# carbon gain at the plant's target nitrogen:carbon ratio, and is partitioned into leaf/root/storage
# organs each step by the ported allocation primitive (nitrogen-conserving):
#
#   d(N)/dt = max(0, NPP)·target_nc_ratio
#   (leaf_N, root_N, storage_N) = allocate_crop_nitrogen(N, leaf_C, root_C, storage_C)
#
# This is a first-order closure of the crop nitrogen loop that keeps the plant N:C near its target as
# biomass grows. The full demand/uptake kinetics (CropNitrogenDemand/CropNitrogenUptake) coupled to
# the soil mineral-N pools, and the Vcmax nitrogen feedback into photosynthesis, are the next
# refinements (the scalar physics for all of these is already ported and tested).

"""
    $(TYPEDEF)

Prognostic crop nitrogen pool and organ partitioning.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropNitrogen{NF} <: Terrarium.AbstractProcess{NF}
    "Organ nitrogen-allocation ratios"
    allocation::CropNitrogenAllocation{NF} = CropNitrogenAllocation(NF)
    "Target plant nitrogen:carbon ratio (gN/gC) governing uptake per unit carbon gain"
    target_nc_ratio::NF = 1 / 30
    "Minimum (structural) leaf nitrogen:carbon ratio — nitrogen limitation reaches 0 here"
    ncleaf_min::NF = 1 / 58.8
    "Reference leaf nitrogen:carbon ratio — nitrogen limitation reaches 1 here"
    ncleaf_ref::NF = 1 / 25
    "Nitrogen turnover rate to litter (per day)"
    turnover_rate::NF = 0.01
end

CropNitrogen(::Type{NF}; kwargs...) where {NF} = CropNitrogen{NF}(; kwargs...)

Terrarium.variables(::CropNitrogen{NF}) where {NF} = (
    Terrarium.prognostic(:crop_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:leaf_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:root_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:storage_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:nitrogen_limitation, XY()),
    Terrarium.auxiliary(:crop_nitrogen_uptake, XY(), units = u"kg/m^2/s"),
    Terrarium.auxiliary(:crop_litterfall_nitrogen, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:net_primary_production, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:leaf_carbon, XY(), units = u"kg/m^2"),
    Terrarium.input(:root_carbon, XY(), units = u"kg/m^2"),
    Terrarium.input(:storage_carbon, XY(), units = u"kg/m^2"),
)

"""
    $(TYPEDSIGNATURES)

Leaf-nitrogen limitation factor ∈ [0,1] on the Rubisco capacity, from the leaf N:C ratio between the
structural minimum (`ncleaf_min` → 0) and the reference (`ncleaf_ref` → 1). Returns 1 when there is no
leaf carbon yet, so early growth is not deadlocked.
"""
@inline function leaf_nitrogen_limitation(n::CropNitrogen{NF}, leaf_nitrogen::NF, leaf_carbon::NF) where {NF}
    nc = leaf_nitrogen / max(leaf_carbon, eps(NF))
    limited = clamp((nc - n.ncleaf_min) / (n.ncleaf_ref - n.ncleaf_min), zero(NF), one(NF))
    return ifelse(leaf_carbon > zero(NF), limited, one(NF))
end

# ---- interface methods --------------------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, n::CropNitrogen, args...)
    out = Terrarium.auxiliary_fields(state, n)
    fields = get_fields(state, n; except = out)
    launch!(grid, XY, compute_crop_nitrogen_auxiliary_kernel!, out, fields, n)
    return nothing
end

""" $(TYPEDSIGNATURES) Acquire nitrogen in proportion to net carbon gain at the target N:C ratio. """
function Terrarium.compute_tendencies!(state, grid, n::CropNitrogen, args...)
    tend = Terrarium.tendency_fields(state, n)
    fields = get_fields(state, n)
    launch!(grid, XY, compute_crop_nitrogen_tendency_kernel!, tend, fields, n)
    return nothing
end

@kernel inbounds = true function compute_crop_nitrogen_auxiliary_kernel!(out, grid, fields, n::CropNitrogen)
    i, j = @index(Global, NTuple)
    leaf, root, storage, _pool = allocate_crop_nitrogen(
        n.allocation, fields.crop_nitrogen[i, j],
        fields.leaf_carbon[i, j], fields.root_carbon[i, j], fields.storage_carbon[i, j], zero(eltype(out.leaf_nitrogen)),
    )
    out.leaf_nitrogen[i, j, 1] = leaf
    out.root_nitrogen[i, j, 1] = root
    out.storage_nitrogen[i, j, 1] = storage
    out.nitrogen_limitation[i, j, 1] = leaf_nitrogen_limitation(n, leaf, fields.leaf_carbon[i, j])
    # Uptake demand (per unit net carbon gain at the target N:C) drawn from the soil mineral N,
    # and turnover of plant nitrogen returned to the soil, both consumed by the soil biogeochemistry.
    NF = eltype(out.crop_nitrogen_uptake)
    out.crop_nitrogen_uptake[i, j, 1] = max(zero(NF), fields.net_primary_production[i, j]) * n.target_nc_ratio
    out.crop_litterfall_nitrogen[i, j, 1] =
        n.turnover_rate / Terrarium.seconds_per_day(NF) * max(zero(NF), fields.crop_nitrogen[i, j])
end

@kernel inbounds = true function compute_crop_nitrogen_tendency_kernel!(tend, grid, fields, n::CropNitrogen)
    i, j = @index(Global, NTuple)
    # d(N)/dt = uptake − litterfall (plant nitrogen: gained by root uptake, lost to soil litter).
    tend.crop_nitrogen[i, j, 1] = fields.crop_nitrogen_uptake[i, j] - fields.crop_litterfall_nitrogen[i, j]
end
