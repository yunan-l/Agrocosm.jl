# Crop nitrogen allocation (LPJmL `crop_nitrogen`). The total plant nitrogen stock is redistributed
# among the leaf, root, storage, and mobile-pool organs in proportion to each organ's carbon divided
# by its target C:N ratio (relative to leaf). This conserves total plant nitrogen. Tested scalar
# physics; assembled in the crop nitrogen coupling (plan Phase 5).

"""
    $(TYPEDEF)

Relative organ C:N ratios for crop nitrogen allocation (leaf is the reference, ratio 1).

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropNitrogenAllocation{NF}
    "Root-to-leaf relative C:N ratio"
    @param ratio_root::NF = 1.16 (bounds = Positive,)
    "Storage-organ-to-leaf relative C:N ratio"
    @param ratio_storage::NF = 0.99 (bounds = Positive,)
    "Mobile-pool-to-leaf relative C:N ratio"
    @param ratio_pool::NF = 3.0 (bounds = Positive,)
end

CropNitrogenAllocation(::Type{NF}; kwargs...) where {NF} = CropNitrogenAllocation{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Allocate the total plant nitrogen among organs by carbon-to-target-C:N weights. Returns
`(leaf_n, root_n, storage_n, pool_n)`, whose sum equals `total_nitrogen` (nitrogen is conserved).
Returns all zeros when there is no organ carbon.
"""
@inline function allocate_crop_nitrogen(
        a::CropNitrogenAllocation{NF}, total_nitrogen::NF,
        leaf_carbon::NF, root_carbon::NF, storage_carbon::NF, pool_carbon::NF,
    ) where {NF}
    leaf_weight = leaf_carbon
    root_weight = root_carbon / a.ratio_root
    storage_weight = storage_carbon / a.ratio_storage
    pool_weight = pool_carbon / a.ratio_pool
    total_weight = leaf_weight + root_weight + storage_weight + pool_weight
    scale = ifelse(total_weight > zero(NF), total_nitrogen / total_weight, zero(NF))
    return leaf_weight * scale, root_weight * scale, storage_weight * scale, pool_weight * scale
end
