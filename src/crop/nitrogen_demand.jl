# Crop nitrogen demand (LPJmL `ndemand_crop`). Leaf demand is the Rubisco nitrogen requirement
# implied by the photosynthetic capacity (temperature-scaled — the inverse of the Vcmax nitrogen
# limitation) plus the structural minimum leaf N; total demand adds the other organs (root, mobile
# pool, storage) at the leaf N:C ratio scaled by their relative C:N ratios. This is the tested
# scalar physics; it is wired into the crop nitrogen coupling.

"""
    $(TYPEDEF)

Parameters of the crop nitrogen demand.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrogenDemand{NF}
    "Atmospheric-pressure/psychrometric scaling constant (LPJmL `p`)"
    @param pressure_scale::NF = 25.0 (bounds = Positive,)
    "Temperature sensitivity of the nitrogen requirement for Rubisco activity"
    @param k_temp::NF = 0.0693 (bounds = Positive,)
    "Minimum leaf nitrogen-to-carbon ratio"
    @param ncleaf_min::NF = 1 / 58.8 (bounds = Positive,)
    "Maximum leaf nitrogen-to-carbon ratio"
    @param ncleaf_max::NF = 1 / 14.3 (bounds = Positive,)
    "Root-to-leaf relative C:N ratio"
    @param ratio_root::NF = 1.16 (bounds = Positive,)
    "Storage-organ-to-leaf relative C:N ratio"
    @param ratio_storage::NF = 0.99 (bounds = Positive,)
    "Mobile-pool-to-leaf relative C:N ratio"
    @param ratio_pool::NF = 3.0 (bounds = Positive,)
end

CropNitrogenDemand(::Type{NF}; kwargs...) where {NF} = CropNitrogenDemand{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Crop leaf and total nitrogen demand from the photosynthetic capacity `vcmax`, the organ carbon
pools, and air temperature (°C). Returns `(nitrogen_demand_leaf, nitrogen_demand_total)`.
"""
@inline function crop_nitrogen_demand(
        d::CropNitrogenDemand{NF}, vcmax::NF,
        leaf_carbon::NF, root_carbon::NF, pool_carbon::NF, storage_carbon::NF, temperature::NF,
    ) where {NF}
    rubisco_demand = d.pressure_scale * NF(1.0e-3) * vcmax /
        (NF(86400) * NF(12) * NF(1.0e-6)) * exp(-d.k_temp * (temperature - NF(25)))
    demand_leaf = rubisco_demand + d.ncleaf_min * leaf_carbon
    nc_ratio = ifelse(leaf_carbon > zero(NF), demand_leaf / leaf_carbon, zero(NF))
    nc_ratio = clamp(nc_ratio, d.ncleaf_min, d.ncleaf_max)
    demand_total = demand_leaf + nc_ratio *
        (root_carbon / d.ratio_root + pool_carbon / d.ratio_pool + storage_carbon / d.ratio_storage)
    return demand_leaf, demand_total
end
