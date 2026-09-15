"""Persistent seasonal crop-water memory."""
mutable struct CropWaterState{A}
    demand_sum::A # Seasonal accumulated transpiration demand (mm).
    supply_sum::A # Seasonal accumulated transpiration supply (mm).
    sufficiency::A # Prior-day water-sufficiency multiplier needed by next-day canopy growth (0–1).
    # Root-weighted soil matric potential, bar and negative. `sufficiency` is the
    # stomatal stress and carries FAO-56's plateau and a clip at one, both right
    # for transpiration and both wrong for expansive growth. What replaced the
    # ratio here is a POTENTIAL and not a content: the same relative content is
    # -1.2 bar in Maricopa's clay loam and -0.4 in Braunschweig's loamy sand, and
    # the crop's own canopy thermometer ranks those two the way the potential
    # does and the content does not.
    root_zone_potential::A
end

"""Current-day crop-water fluxes."""
mutable struct CropWaterFluxes{A, M}
    interception::A       # Rainfall intercepted and evaporated by the canopy (mm day⁻¹).
    transpiration_layer::M # Root-water uptake/transpiration from each soil layer (mm day⁻¹).
end

function init_crop_water_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropWaterState(float_state(), float_state(), device(ones(T, cell_size)),
                          device(fill(T(-1) / T(3), cell_size)))
end

function init_crop_water_fluxes(::Type{T}, cell_size::Int, device;
                                soil_layers::Int = 5) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropWaterFluxes(
        float_flux(),
        device(zeros(T, soil_layers, cell_size)),
    )
end
