"""Daily crop photosynthesis fluxes, capacities, and stress state."""
mutable struct CropPhotosynthesis{A}
    gross_assimilation::A
    net_assimilation::A
    water_limited_assimilation::A
    leaf_respiration::A
    vmax::A
    lambda::A
    temperature_stress::A
end

function init_crop_photosynthesis(cell_size::Int, device)
    float_state() = device(zeros(Float32, cell_size))
    return CropPhotosynthesis(ntuple(_ -> float_state(), 7)...)
end
