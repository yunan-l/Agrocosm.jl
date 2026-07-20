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
function init_crop(cell_size::Int,
                   device;
                   carbon_pools::Int = 4,
                   soil_layers::Int = 5)

    crop = Crop(
        init_crop_phenology(cell_size, device),
        init_crop_canopy(cell_size, device),
        init_crop_carbon(cell_size, device; carbon_pools = carbon_pools),
        init_crop_nitrogen(cell_size, device),
        init_crop_water(cell_size, device; soil_layers = soil_layers),
        init_crop_calendar(cell_size, device),
        init_crop_photosynthesis(cell_size, device),
    )

    return crop
end
