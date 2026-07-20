"""Crop water fluxes, root distribution, and water-stress state."""
mutable struct CropWater{A, M}
    canopy_conductance::A
    transpiration::A
    canopy_wet::A
    interception::A
    transpiration_layer::M
    root_distribution::A
    deficit::A
    demand_sum::A
    supply_sum::A
    stress::A
    waterlogging_days::A
    waterlogging_stress::A
    root_zone_water::A
end

init_crop_water(cell_size::Int, device; kwargs...) =
    init_crop_water(Float32, cell_size, device; kwargs...)
function init_crop_water(::Type{T}, cell_size::Int, device;
                         soil_layers::Int = 5) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))

    return CropWater(
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        device(zeros(T, soil_layers, cell_size)),
        device(zeros(T, soil_layers)),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
        float_state(),
    )
end
