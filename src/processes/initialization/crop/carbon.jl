"""Persistent plant carbon stocks."""
mutable struct CropCarbonState{A}
    biomass::A # Total living crop carbon biomass (gC m⁻²).
    leaf::A    # Living leaf carbon stock (gC m⁻²).
    root::A    # Living root carbon stock (gC m⁻²).
    pool::A    # Mobile/intermediate plant carbon pool (gC m⁻²).
    storage::A # Harvestable storage-organ carbon stock (gC m⁻²).
end

"""Current-day plant carbon fluxes."""
mutable struct CropCarbonFluxes{A}
    yield::A                      # Grain/storage carbon harvested today (gC m⁻² day⁻¹).
    harvest_export::A             # Total crop carbon removed at harvest (gC m⁻² day⁻¹).
    npp::A                        # Net primary production after all plant respiration (gC m⁻² day⁻¹).
    respiration::A                # Root, organ-maintenance, and growth respiration (gC m⁻² day⁻¹).
    gross_assimilation::A         # Gross daily canopy carbon assimilation/GPP (gC m⁻² day⁻¹).
    net_assimilation::A           # Nonnegative daytime assimilation after leaf respiration (gC m⁻² day⁻¹).
    water_limited_assimilation::A # Assimilation expressed as transpiration demand (mm day⁻¹).
    leaf_respiration::A           # Daily leaf/dark respiration (gC m⁻² day⁻¹).
    biological_fixation_cost::A   # Carbon allocated to biological N fixation (gC m⁻² day⁻¹).
end

function init_crop_carbon_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropCarbonState(ntuple(_ -> float_state(), 5)...)
end

function init_crop_carbon_fluxes(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropCarbonFluxes(ntuple(_ -> float_flux(), 9)...)
end
