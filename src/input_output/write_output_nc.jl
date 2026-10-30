"""
write_output_nc(output_path, varname, values; ...)

Write model outputs to NetCDF on the longitude/latitude grid defined by `grid_path`.
Provide one `(longitude_index, latitude_index)` pair in `grid_indices` for each
input cell. When `time` is omitted for a time series, use `1:ntime`.

`values` can be:
- `AbstractVector` of length `ncell` (single snapshot), or
- `AbstractMatrix` of size `(ntime, ncell)` (time series).

`grid_indices = [(lon_idx, lat_idx), ...]`.
"""
function write_output_nc(
    output_path::String,
    varname::String,
    values::AbstractArray;
    grid_path::String,
    grid_indices::AbstractVector{<:Tuple{<:Integer, <:Integer}},
    time::Union{Nothing, AbstractVector} = nothing,
    units::String = "",
    long_name::String = "",
    fill_value::Float32 = -9999.0f0,
)
    ndims(values) in (1, 2) || throw(ArgumentError(
        "values must be 1D (cell) or 2D (time, cell)",
    ))

    vals = Array(values)
    lat, lon, cell_map = _read_grid_spec(grid_path)
    ntime = ndims(vals) == 1 ? 1 : size(vals, 1)
    ncell = ndims(vals) == 1 ? length(vals) : size(vals, 2)
    length(grid_indices) == ncell || throw(DimensionMismatch(
        "grid_indices has $(length(grid_indices)) entries for $ncell input cells",
    ))
    all(index -> 1 <= index[1] <= length(lon) && 1 <= index[2] <= length(lat), grid_indices) ||
        throw(ArgumentError("grid_indices contains a coordinate outside the grid"))
    time !== nothing && length(time) != ntime && throw(DimensionMismatch(
        "time has $(length(time)) entries for $ntime output rows",
    ))
    cell_map = ifelse.(ismissing.(cell_map), missing, 0)

    if ntime > 1
        init_out = repeat(cell_map; outer=(1, 1, ntime))
    else
        init_out = cell_map
    end

    isfile(output_path) && rm(output_path; force=true)
    NCDataset(output_path, "c") do ds
        defDim(ds, "longitude", length(lon))
        defDim(ds, "latitude", length(lat))
        if ntime > 1 || time !== nothing
            defDim(ds, "time", ntime)
        end

        vlon = defVar(ds, "longitude", Float32, ("longitude",))
        vlat = defVar(ds, "latitude", Float32, ("latitude",))
        vlon[:] = Float32.(lon)
        vlat[:] = Float32.(lat)
        vlon.attrib["units"] = "degrees_east"
        vlat.attrib["units"] = "degrees_north"

        if ntime > 1 || time !== nothing
            vtime = defVar(ds, "time", Float32, ("time",))
            vtime[:] = Float32.(time === nothing ? (1:ntime) : time)
            vtime.attrib["long_name"] = "time index"
        end

        dims = (ntime > 1 || time !== nothing) ? ("longitude", "latitude", "time") : ("longitude", "latitude")
        vout = defVar(ds, varname, Float32, dims; fillvalue = fill_value)
        
        if !isempty(units)
            vout.attrib["units"] = units
        end
        if !isempty(long_name)
            vout.attrib["long_name"] = long_name
        end

        if ntime == 1
            vout[:, :] .= init_out
            for (i, grid_idx) in enumerate(grid_indices)
                lonidx, latidx = grid_idx
                vout[lonidx, latidx] .= vals[i]
            end
        else
            vout[:, :, :] .= init_out
            for (i, grid_idx) in enumerate(grid_indices)
                lonidx, latidx = grid_idx
                vout[lonidx, latidx, :] .= vals[:, i]
            end
        end
    end

    return output_path
end

function _read_grid_spec(grid_path::String)
    ds = NCDataset(grid_path, "r")
    lon = ds["longitude"][:]
    lat = ds["latitude"][:]
    cell_map = ds["cellid"][:, :]
    close(ds)
    return lat, lon, cell_map
end
