"""Weather controls, in model units, along the third forcing-array dimension."""
const WEATHER_VARIABLES = (:temp, :prec, :sw, :lw, :wind)

"""Pack prepared daily weather as `(day, cell, variable)`; CO₂ and N deposition stay fixed."""
function weather_forcing(climate::NamedTuple)
    arrays = map(name -> getproperty(climate, name), WEATHER_VARIABLES)
    all(array -> ndims(array) == 2, arrays) ||
        throw(DimensionMismatch("weather fields must be (day, cell) matrices"))
    all(array -> size(array) == size(first(arrays)), arrays) ||
        throw(DimensionMismatch("weather fields must have identical (day, cell) shapes"))
    return cat(arrays...; dims = 3)
end

# A distinct opt-in input boundary: the production readclimate! rule remains
# inactive for existing parameter objectives. Both CPU and GPU use this kernel.
function apply_weather_forcing!(weather::DailyWeather, forcing, day::Integer)
    ndims(forcing) == 3 && size(forcing, 3) == length(WEATHER_VARIABLES) &&
        size(forcing, 2) == length(weather.temp) || throw(DimensionMismatch("weather control shape mismatch"))
    1 <= day <= size(forcing, 1) || throw(BoundsError(forcing, (day, :, :)))
    launch_1D!(
        apply_weather_forcing_kernel!, weather.temp, weather.prec, weather.swr,
        weather.lwr, weather.wind, forcing, day,
    )
    return nothing
end

@kernel inbounds = true function apply_weather_forcing_kernel!(
    temperature, precipitation, shortwave, longwave, wind, forcing, day,
)
    cell = @index(Global)
    temperature[cell] = forcing[day, cell, 1]
    precipitation[cell] = forcing[day, cell, 2]
    shortwave[cell] = forcing[day, cell, 3]
    longwave[cell] = forcing[day, cell, 4]
    wind[cell] = forcing[day, cell, 5]
end

"""Ordinary one-cell production replay through a specified harvest; requires Enzyme extension."""
function weather_harvest_replay end

"""Fixed-harvest weather-to-yield reverse gradient; requires `using Enzyme`."""
function enzyme_weather_harvest_gradient end

"""Fixed-harvest weather directional derivative; requires `using Enzyme`."""
function enzyme_weather_forward_directional end

"""Fixed-harvest reverse gradient with respect to named `CFTParameters` fields.

The process-attribution counterpart of `enzyme_weather_harvest_gradient`:
`d(yield)/d(theta)` for the parameters that define each stress mechanism.
Requires `using Enzyme`.
"""
function enzyme_process_parameter_gradient end
