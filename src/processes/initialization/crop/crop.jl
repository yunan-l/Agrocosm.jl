"""Process-grouped crop state; leaf fields remain backend arrays (CPU/GPU SoA)."""
mutable struct Crop{P, C, B, N, W, K, F}
    phenology::P
    canopy::C
    carbon::B
    nitrogen::N
    water::W
    calendar::K
    photosynthesis::F
end

"""
init_crop(cell_size, device; carbon_pools=4, soil_layers=5)

Allocate the complete process-grouped crop state on `device`. Calendar and
photosynthesis state are owned by `Crop`.
"""
init_crop(cell_size::Int, device; kwargs...) =
    init_crop(Float32, cell_size, device; kwargs...)
function init_crop(::Type{T},
                   cell_size::Int,
                   device;
                   carbon_pools::Int = 4,
                   soil_layers::Int = 5) where {T <: AbstractFloat}

    crop = Crop(
        init_crop_phenology(T, cell_size, device),
        init_crop_canopy(T, cell_size, device),
        init_crop_carbon(T, cell_size, device; carbon_pools = carbon_pools),
        init_crop_nitrogen(T, cell_size, device),
        init_crop_water(T, cell_size, device; soil_layers = soil_layers),
        init_crop_calendar(cell_size, device),
        init_crop_photosynthesis(T, cell_size, device),
    )

    return crop
end
