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
    "Leaf-nitrogen Rubisco-capacity limitation (the LPJmL Vcmax nitrogen feedback into photosynthesis)"
    vcmax_limit::CropNitrogenVcmaxLimit{NF} = CropNitrogenVcmaxLimit(NF)
    "Crop nitrogen demand (from photosynthetic capacity + organ pools)"
    demand::CropNitrogenDemand{NF} = CropNitrogenDemand(NF)
    "Root nitrogen uptake kinetics (Michaelis–Menten from the soil mineral-N pools)"
    uptake::CropNitrogenUptakeKinetics{NF} = CropNitrogenUptakeKinetics(NF)
end

CropNitrogen(::Type{NF}; kwargs...) where {NF} = CropNitrogen{NF}(; kwargs...)

Terrarium.variables(::CropNitrogen{NF}) where {NF} = (
    Terrarium.prognostic(:crop_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:leaf_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:root_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:storage_nitrogen, XY(), units = u"kg/m^2"),
    Terrarium.auxiliary(:nitrogen_capacity, XY(), units = u"g/m^2/s"),
    Terrarium.auxiliary(:crop_nitrogen_uptake, XY(), units = u"kg/m^2/s"),
    Terrarium.auxiliary(:crop_litterfall_nitrogen, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:net_primary_production, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:leaf_carbon, XY(), units = u"kg/m^2"),
    Terrarium.input(:root_carbon, XY(), units = u"kg/m^2"),
    Terrarium.input(:storage_carbon, XY(), units = u"kg/m^2"),
    Terrarium.input(:air_temperature, XY(), default = NF(10), units = u"°C"),
    # Demand + uptake kinetics couplings: the potential Vcmax (photosynthesis) sizes the demand, and the
    # soil mineral-N pools + root fraction + soil temperature (over the root zone) limit the supply.
    Terrarium.input(:potential_vcmax, XY(), default = zero(NF), units = u"g/m^2/s"),
    Terrarium.input(:soil_ammonium, XYZ(), default = zero(NF), units = u"kg/m^3"),
    Terrarium.input(:soil_nitrate, XYZ(), default = zero(NF), units = u"kg/m^3"),
    Terrarium.input(:root_fraction, XYZ(), default = zero(NF)),
    Terrarium.input(:temperature, XYZ(), default = NF(5), units = u"°C"),
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
    NF = eltype(out.crop_nitrogen_uptake)
    seconds_per_day = Terrarium.seconds_per_day(NF)
    leaf_carbon = fields.leaf_carbon[i, j]
    root_carbon = fields.root_carbon[i, j]
    storage_carbon = fields.storage_carbon[i, j]
    leaf, root, storage, _pool = allocate_crop_nitrogen(
        n.allocation, fields.crop_nitrogen[i, j], leaf_carbon, root_carbon, storage_carbon, zero(NF),
    )
    out.leaf_nitrogen[i, j, 1] = leaf
    out.root_nitrogen[i, j, 1] = root
    out.storage_nitrogen[i, j, 1] = storage
    # Nitrogen-supported Rubisco capacity (gC/m²/s): photosynthesis caps its light-derived Vc_max at
    # this, applying the LPJmL leaf-nitrogen limitation. Large when there is no leaf carbon yet.
    out.nitrogen_capacity[i, j, 1] = nitrogen_supported_vcmax(n.vcmax_limit, leaf, leaf_carbon, fields.air_temperature[i, j])

    # Nitrogen demand (gN/m²) sized from the potential (light-derived) Vc_max and the organ carbon
    # pools; convert the revised units (Vc_max /s → /day, carbon kg → g).
    _demand_leaf, demand_total = crop_nitrogen_demand(
        n.demand, fields.potential_vcmax[i, j] * seconds_per_day,
        leaf_carbon * NF(1000), root_carbon * NF(1000), zero(NF), storage_carbon * NF(1000), fields.air_temperature[i, j],
    )
    demand_nitrogen = demand_total / NF(1000)   # kgN/m²

    # Root-zone soil mineral nitrogen (kgN/m²) and the root-weighted soil temperature, integrated over
    # the column (the root fraction sums to unity).
    field_grid = get_field_grid(grid)
    available_ammonium = zero(NF)
    available_nitrate = zero(NF)
    root_zone_temperature = zero(NF)
    for k in 1:size(fields.soil_ammonium, 3)
        Δz = Δzᵃᵃᶜ(i, j, k, field_grid)
        available_ammonium += max(zero(NF), fields.soil_ammonium[i, j, k]) * Δz
        available_nitrate += max(zero(NF), fields.soil_nitrate[i, j, k]) * Δz
        root_zone_temperature += fields.temperature[i, j, k] * fields.root_fraction[i, j, k]
    end

    # Michaelis–Menten root uptake potential per mineral pool (gN/m²/day); the root factor is the
    # soil-temperature response × root carbon. Convert to a per-second rate.
    root_factor = nitrogen_uptake_temperature_response(root_zone_temperature, NF(-25), NF(15), NF(15)) * root_carbon * NF(1000)
    potential_ammonium = root_nitrogen_uptake_potential(n.uptake, available_ammonium * NF(1000), one(NF), root_factor)
    potential_nitrate = root_nitrogen_uptake_potential(n.uptake, available_nitrate * NF(1000), one(NF), root_factor)
    potential_rate = (potential_ammonium + potential_nitrate) / NF(1000) / seconds_per_day   # kgN/m²/s

    # Uptake fills the demand deficit over ~a day, limited by the soil-supply kinetics.
    demand_rate = max(zero(NF), demand_nitrogen - fields.crop_nitrogen[i, j]) / seconds_per_day
    out.crop_nitrogen_uptake[i, j, 1] = min(demand_rate, potential_rate)

    # Turnover of plant nitrogen returned to the soil litter.
    out.crop_litterfall_nitrogen[i, j, 1] =
        n.turnover_rate / seconds_per_day * max(zero(NF), fields.crop_nitrogen[i, j])
end

@kernel inbounds = true function compute_crop_nitrogen_tendency_kernel!(tend, grid, fields, n::CropNitrogen)
    i, j = @index(Global, NTuple)
    # d(N)/dt = uptake − litterfall (plant nitrogen: gained by root uptake, lost to soil litter).
    tend.crop_nitrogen[i, j, 1] = fields.crop_nitrogen_uptake[i, j] - fields.crop_litterfall_nitrogen[i, j]
end
