_test_backend(reference) = values -> begin
    destination = similar(reference, eltype(values), size(values))
    copyto!(destination, values)
    return destination
end

"""Build a canonical lifecycle state over existing crop and soil test arrays."""
function test_model_state(
    crop::Agrocosm.Crop,
    soil::Agrocosm.Soil;
    managed_land = nothing,
    pet = nothing,
    climbuf = nothing,
    weather = nothing,
    output = nothing,
)
    cells = length(crop.state.canopy.lai)
    T = eltype(crop.state.canopy.lai)
    reference = crop.state.canopy.lai
    backend = _test_backend(reference)
    managed_land = isnothing(managed_land) ? Agrocosm.init_managed_land(T, cells, backend) : managed_land
    pet = isnothing(pet) ? Agrocosm.init_pet(T, cells, backend) : pet
    climbuf = isnothing(climbuf) ? Agrocosm.init_climbuf(T, cells, backend) : climbuf
    weather = isnothing(weather) ? Agrocosm.init_weather(T, cells, backend) : weather
    output = isnothing(output) ? Agrocosm.init_output(T, cells, backend) : output
    return Agrocosm.model_state(climbuf, crop, pet, soil, managed_land, weather, output)
end

function test_model_state(crop::Agrocosm.Crop; kwargs...)
    T = eltype(crop.state.canopy.lai)
    cells = length(crop.state.canopy.lai)
    backend = _test_backend(crop.state.canopy.lai)
    soil = Agrocosm.init_soil(T, cells, T.(Agrocosm.soilparams.soildepth), backend)
    return test_model_state(crop, soil; kwargs...)
end

function test_model_state(soil::Agrocosm.Soil; kwargs...)
    T = eltype(soil.water.storage)
    cells = size(soil.water.storage, 2)
    backend = _test_backend(soil.water.storage)
    crop = Agrocosm.init_crop(T, cells, backend)
    return test_model_state(crop, soil; kwargs...)
end
