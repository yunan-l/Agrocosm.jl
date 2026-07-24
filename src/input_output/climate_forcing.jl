# Building Terrarium time-varying input sources from tabulated climate series. This reuses
# Oceananigans `FieldTimeSeries` and Terrarium's `FieldTimeSeriesInputSource`: each series is sampled
# at the given `times` and interpolated to the current clock time on every step (Terrarium's
# `update_inputs!`, invoked from `update_state!` inside `run!`/`timestep!`). It is format-agnostic —
# the caller supplies plain vectors, whether read from JLD2, NetCDF, or generated.

"""
    $(TYPEDSIGNATURES)

Build a Terrarium `InputSources` collection of time-varying **surface** climate inputs from tabulated
series. Each keyword `name = data` (with `data` a vector the same length as `times`) becomes a surface
(`XY`) Oceananigans `FieldTimeSeries` sampled at `times` (seconds since the simulation start), wrapped
in a `FieldTimeSeriesInputSource`. Passing the result to `initialize(model; inputs = …)` drives the
model with the forcing: Terrarium interpolates each series to the current clock time on every step.

The series are horizontally uniform (one value per time, broadcast over the grid columns), covering the
single-column and shared-forcing cases. Each `name` must match a model input variable, e.g.
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
        length(data) == length(times) ||
            throw(ArgumentError("series `$name` has $(length(data)) values but `times` has $(length(times))"))
        fts = FieldTimeSeries(grid, XY(), times)
        snapshots = interior(fts)                     # (Nx, Ny, Nz, Nt)
        for n in eachindex(times)
            @inbounds snapshots[:, :, :, n] .= data[n]
        end
        unit = get(all_units, name, nothing)
        return isnothing(unit) ? InputSource(fts; name = name) : InputSource(fts; name = name, units = unit)
    end
    return InputSources(sources...)
end
