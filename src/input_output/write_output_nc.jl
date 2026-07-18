"""
write_output_nc(output_path, varname, values; ...)

Write model outputs to NetCDF on the longitude/latitude grid defined by `grid_path`
(default: `inputs/grid.nc`).

`values` can be:
- `AbstractVector` of length `ncell` (single snapshot), or
- `AbstractMatrix` of size `(ntime, ncell)` (time series).

By default, columns are interpreted as global LPJmL-style cell ids in ascending order.
`grid_indices = [(lon_idx, lat_idx), ...]`.
"""
function write_output_nc(
    output_path::String,
    varname::String,
    values::AbstractArray;
    grid_path::String = joinpath(pkgdir(Agrocosm), "inputs/grid.nc"),
    grid_indices::Union{Nothing, AbstractVector{<:Tuple{<:Integer, <:Integer}}} = nothing,
    time::Union{Nothing, AbstractVector} = nothing,
    units::String = "",
    long_name::String = "",
    fill_value::Float32 = -9999.0f0,
)
    @assert ndims(values) in (1, 2) "values must be 1D (cell) or 2D (time, cell)"

    vals = Array(values)
    lat, lon, cell_map = _read_grid_spec(grid_path)
    cell_map = ifelse.(ismissing.(cell_map), missing, 0)

    ntime = ndims(vals) == 1 ? 1 : size(vals, 1)
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
            vtime[:] = Float32.(time)
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