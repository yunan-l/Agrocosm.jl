"""Persistent crop state required to advance the next daily transition."""
mutable struct CropState{P, C, B, N, W, K}
    phenology::P
    canopy::C
    carbon::B
    nitrogen::N
    water::W
    calendar::K
end

"""Current-day carbon, nitrogen, and water transfers."""
mutable struct CropFluxes{C, N, W}
    carbon::C
    nitrogen::N
    water::W
end

"""Recomputable current-day canopy, photosynthesis, and stress variables."""
mutable struct CropStressAuxiliary{A, M}
    nitrogen_demand_total::A
    nitrogen_demand_leaf::A
    nitrogen::A
    nitrogen_deficit::A
    water_deficit::A
    water::A
    waterlogging::A
    root_zone_water::A
    root_distribution::M
end

mutable struct CropAuxiliary{C, P, S}
    canopy::C
    photosynthesis::P
    stress::S
end

"""Preallocated implementation-only buffers, excluded from restart and output."""
mutable struct CropWorkspace{A}
    respiration_temperature_response::A
end

"""Lifecycle-grouped crop storage; all leaves remain backend arrays (CPU/GPU SoA)."""
mutable struct Crop{S, F, A, E, W}
    state::S
    fluxes::F
    auxiliary::A
    events::E
    workspace::W
end

"""
init_crop(cell_size, device; soil_layers=5)

Allocate the complete lifecycle-grouped crop storage on `device`. Calendar
state and photosynthesis auxiliary variables are owned by `Crop`.
"""
init_crop(cell_size::Int, device; kwargs...) =
    init_crop(Float32, cell_size, device; kwargs...)
function init_crop(::Type{T},
                   cell_size::Int,
                   device;
                   soil_layers::Int = 5) where {T <: AbstractFloat}

    float_auxiliary() = device(zeros(T, cell_size))
    crop = Crop(
        CropState(
            init_crop_phenology(T, cell_size, device),
            init_crop_canopy_state(T, cell_size, device),
            init_crop_carbon_state(T, cell_size, device),
            init_crop_nitrogen_state(T, cell_size, device),
            init_crop_water_state(T, cell_size, device),
            init_crop_calendar_state(cell_size, device),
        ),
        CropFluxes(
            init_crop_carbon_fluxes(T, cell_size, device),
            init_crop_nitrogen_fluxes(T, cell_size, device),
            init_crop_water_fluxes(T, cell_size, device; soil_layers = soil_layers),
        ),
        CropAuxiliary(
            init_crop_canopy_auxiliary(T, cell_size, device),
            init_crop_photosynthesis_auxiliary(T, cell_size, device),
            CropStressAuxiliary(
                float_auxiliary(), float_auxiliary(), float_auxiliary(),
                float_auxiliary(), float_auxiliary(), float_auxiliary(),
                float_auxiliary(), float_auxiliary(), device(zeros(T, soil_layers)),
            ),
        ),
        init_crop_events(cell_size, device),
        CropWorkspace(float_auxiliary()),
    )

    return crop
end
