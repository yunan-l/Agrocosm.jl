"""Persistent plant carbon stocks."""
mutable struct CropCarbonState{A}
    biomass::A
    leaf::A
    root::A
    pool::A
    storage::A
end

"""Current-day plant carbon fluxes."""
mutable struct CropCarbonFluxes{A}
    yield::A
    npp::A
    respiration::A
    gross_assimilation::A
    net_assimilation::A
    water_limited_assimilation::A
    leaf_respiration::A
end

function init_crop_carbon_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropCarbonState(ntuple(_ -> float_state(), 5)...)
end

function init_crop_carbon_fluxes(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropCarbonFluxes(ntuple(_ -> float_flux(), 7)...)
end
