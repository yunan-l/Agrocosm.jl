"""Layered soil water stocks, hydraulic properties, and daily fluxes."""
mutable struct SoilWater{A, M}
    storage::M                  # Liquid soil-water stock in each layer (mm).
    ice_storage::M              # Total layer ice stock; sum of the three LPJmL ice pools (mm).
    wilting_ice_fraction::M     # Fraction of permanent-wilting water frozen as `ice_pwp` (0–1).
    available_ice_storage::M    # Ice within the plant-available interval, `ice_depth` (mm).
    free_ice_storage::M         # Ice within gravitational/free water, `ice_fw` (mm).
    evaporation::M              # Soil evaporation removed from each layer today (mm day⁻¹).
    relative_content::M         # Relative plant-available liquid water content (0–1).
    free_water::M               # Liquid water above field capacity (mm).
    wilting_fraction::M         # Permanent-wilting water as a fraction of layer depth (0–1).
    wilting_storage::M          # Water storage at permanent wilting point (mm).
    field_capacity::M           # Volumetric/relative field-capacity threshold (fraction).
    saturation_fraction::M      # Saturated volumetric water content (fraction).
    saturation_storage::M       # Total water storage at saturation (mm).
    beta::M                     # Clapp–Hornberger soil-water retention exponent (dimensionless).
    holding_capacity_fraction::M # Plant-available holding capacity as fraction of depth (0–1).
    holding_capacity_storage::M # Plant-available water-holding capacity (mm).
    saturated_conductivity::M   # Saturated hydraulic conductivity (mm day⁻¹).
    influx::M                   # Water entering each layer during current infiltration step (mm day⁻¹).
    outflux::M                  # Water leaving each layer during current infiltration step (mm day⁻¹).
    surface_runoff::A           # Surface runoff leaving the soil column (mm day⁻¹).
    lateral_runoff::M           # Lateral drainage leaving each soil layer (mm day⁻¹).
    bottom_drainage::A          # Drainage leaving the bottom soil layer (mm day⁻¹).
    infiltration::A             # Rain/melt water remaining for soil infiltration (mm day⁻¹).
    percolation::M              # Downward layer-to-layer percolation (mm day⁻¹).
    # Water the surface refused, held instead of deleted. LPJmL removes rejected
    # rainfall as surface runoff on the day it arrives, which is why this lineage
    # cannot waterlog at all: measured on five real cells, root-zone water-filled
    # pore space never reaches 0.90, not even irrigated paddy rice on 54% clay,
    # so the aeration-stress band the rest of the field uses (0.85-0.95) is
    # unreachable. See `docs/32`. This store is offered back to infiltration the
    # next day, so it is a delay rather than a new source, and it is capped by
    # `ponding_capacity`, which ships at zero and reproduces LPJmL bitwise.
    ponding::A                  # Water ponded on the surface, carried between days (mm).
end

init_soil_water(cell_size::Int, device; kwargs...) =
    init_soil_water(Float32, cell_size, device; kwargs...)
function init_soil_water(::Type{T}, cell_size::Int, device;
                         soil_layers::Int = 5) where {T <: AbstractFloat}
    layer_state() = device(zeros(T, soil_layers, cell_size))
    cell_state() = device(zeros(T, cell_size))
    return SoilWater(
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        cell_state(),
        layer_state(),
        cell_state(),
        cell_state(),
        layer_state(),
        cell_state(),
    )
end
