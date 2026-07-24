# Building Terrarium time-varying input sources from tabulated climate series. This reuses
# Oceananigans `FieldTimeSeries` and Terrarium's `FieldTimeSeriesInputSource`: each series is sampled
# at the given `times` and interpolated to the current clock time on every step (Terrarium's
# `update_inputs!`, invoked from `update_state!` inside `run!`/`timestep!`). It is format-agnostic —
# the caller supplies plain vectors, whether read from JLD2, NetCDF, or generated.

# Populate a surface `FieldTimeSeries`' snapshots from either a horizontally-uniform series (one value
# per time) or per-column series (a `time × column` matrix).
function _fill_snapshots!(fts, data::AbstractVector)
    snapshots = interior(fts)                         # (Nx, Ny, Nz, Nt)
    for n in axes(snapshots, 4)
        @inbounds snapshots[:, :, :, n] .= data[n]
    end
    return fts
end

function _fill_snapshots!(fts, data::AbstractMatrix)
    snapshots = interior(fts)                         # (Ncolumns, 1, 1, Nt)
    ncolumns, _, _, ntimes = size(snapshots)
    size(data) == (ntimes, ncolumns) ||
        throw(ArgumentError("per-column series must be (ntimes=$ntimes, ncolumns=$ncolumns); got $(size(data))"))
    for n in 1:ntimes, c in 1:ncolumns
        @inbounds snapshots[c, 1, 1, n] = data[n, c]
    end
    return fts
end

series_ntimes(data::AbstractVector) = length(data)
series_ntimes(data::AbstractMatrix) = size(data, 1)

"""
    $(TYPEDSIGNATURES)

Build a Terrarium `InputSources` collection of time-varying **surface** climate inputs from tabulated
series. Each keyword `name = data` becomes a surface (`XY`) Oceananigans `FieldTimeSeries` sampled at
`times` (seconds since the simulation start), wrapped in a `FieldTimeSeriesInputSource`. Passing the
result to `initialize(model; inputs = …)` drives the model with the forcing: Terrarium interpolates each
series to the current clock time on every step.

`data` is either a **vector** the same length as `times` (horizontally uniform — broadcast over all grid
columns, e.g. a single column or shared forcing), or a **`time × column` matrix** giving each column its
own series (for multi-column / global grids). Each `name` must match a model input variable, e.g.
`air_temperature`, `surface_shortwave_down`.

```julia
times = (0:ndays - 1) .* Terrarium.seconds_per_day(Float64)
inputs = surface_climate_inputs(grid, times; air_temperature = temp, surface_shortwave_down = swdown)
integrator = initialize(model; inputs)
```
"""
# Units of the standard surface climate inputs, so the input-source variable declaration matches the
# model's (mismatched units are a duplicate-variable conflict at `initialize`). Override via `units`.
const SURFACE_CLIMATE_UNITS = (
    air_temperature = u"°C",
    surface_shortwave_down = u"W/m^2",
    surface_longwave_down = u"W/m^2",
    air_pressure = u"Pa",
    windspeed = u"m/s",
)

function surface_climate_inputs(grid::Terrarium.AbstractLandGrid, times::AbstractVector; units = (;), series...)
    isempty(series) && throw(ArgumentError("provide at least one named climate series"))
    all_units = merge(SURFACE_CLIMATE_UNITS, units)
    sources = map(collect(pairs(series))) do (name, data)
        series_ntimes(data) == length(times) ||
            throw(ArgumentError("series `$name` has $(series_ntimes(data)) time samples but `times` has $(length(times))"))
        fts = FieldTimeSeries(grid, XY(), times)
        _fill_snapshots!(fts, data)
        unit = get(all_units, name, nothing)
        return isnothing(unit) ? InputSource(fts; name = name) : InputSource(fts; name = name, units = unit)
    end
    return InputSources(sources...)
end
