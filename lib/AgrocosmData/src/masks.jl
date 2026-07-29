function _valid_crop_fraction(value, threshold)
    ismissing(value) && return false
    value isa AbstractFloat && !isfinite(value) && throw(ArgumentError("landuse contains a non-finite value"))
    return value > threshold
end

"""
Build a fixed allocation selection from the temporal union of positive land use.
The returned `active` matrix retains annual activity only for allocated cells.
"""
function build_crop_mask(
    grid::GridIndex,
    landuse::AbstractMatrix;
    selection::CellSelection = all_cells(grid),
    threshold::Real = 0,
)
    size(landuse, 2) == length(selection.cell_ids) ||
        throw(DimensionMismatch("landuse must have dimensions time × selected cell"))
    active_full = BitMatrix(undef, size(landuse))
    for index in eachindex(landuse)
        active_full[index] = _valid_crop_fraction(landuse[index], threshold)
    end
    allocated_local = findall(vec(any(active_full; dims = 1)))
    isempty(allocated_local) && throw(ArgumentError("selected PFT has no active land-use cells"))
    compact_indices = selection.compact_indices[allocated_local]
    allocated = select_cells(grid, compact_indices)

    value_type = Base.nonmissingtype(eltype(landuse))
    fraction = Matrix{value_type}(undef, size(landuse, 1), length(allocated_local))
    for (output_cell, input_cell) in pairs(allocated_local), time_index in axes(landuse, 1)
        value = landuse[time_index, input_cell]
        fraction[time_index, output_cell] = ismissing(value) ? zero(value_type) : value_type(value)
    end
    return CropMask(allocated, fraction, active_full[:, allocated_local])
end
