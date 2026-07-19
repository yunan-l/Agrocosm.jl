"""Layered soil temperature and thermal-diffusivity parameters."""
mutable struct SoilThermal{A, M}
    temperature::M
    diffusivity_0::A
    diffusivity_15::A
end

function init_soil_thermal(cell_size::Int, device; soil_layers::Int = 5)
    return SoilThermal(
        device(zeros(Float32, soil_layers, cell_size)),
        device(zeros(Float32, cell_size)),
        device(zeros(Float32, cell_size)),
    )
end
