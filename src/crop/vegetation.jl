# A minimal crop vegetation model wiring the crop processes into Terrarium's `AbstractVegetation`
# interface, so it slots into a `LandModel` alongside Terrarium's soil, surface, and atmosphere. The
# crop leaf area index is driven by the prognostic phenological heat units (not a carbon-pool
# equilibrium): heat units accumulate with temperature, set the heat-unit fraction, which drives the
# LAI trajectory, which feeds photosynthesis. It assembles the ported crop primitives with the
# multi-organ carbon/nitrogen pools and the soil C–N coupling.

"""
    $(TYPEDEF)

Minimal crop vegetation model: prognostic phenological heat units → heat-unit-driven LAI →
stomatal conductance → C3/C4 photosynthesis. Wires into a Terrarium `LandModel` as the vegetation
component.

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropVegetation{
        NF,
        PhenologyDynamics <: CropPhenologyDynamics{NF},
        Phenology <: CropPhenology{NF},
        StomatalConductance <: CropStomatalConductance{NF},
        Photosynthesis <: CropPhotosynthesis{NF},
        RootDistribution <: Terrarium.AbstractRootDistribution{NF},
        PlantAvailableWater <: Union{Nothing, Terrarium.AbstractPlantAvailableWater{NF}},
        Carbon <: CropCarbon{NF},
        Nitrogen <: CropNitrogen{NF},
    } <: Terrarium.AbstractVegetation{NF}
    "Phenological heat-unit accumulation (prognostic)"
    phenology_dynamics::PhenologyDynamics
    "Heat-unit-driven leaf-area-index trajectory"
    phenology::Phenology
    "Crop stomatal conductance (λ, canopy conductance)"
    stomatal_conductance::StomatalConductance
    "C3/C4 crop photosynthesis"
    photosynthesis::Photosynthesis
    "Root distribution (root fraction per soil layer)"
    root_distribution::RootDistribution
    "Plant-available-water soil-moisture stress factor β, or `nothing` for a well-watered crop (β=1)"
    plant_available_water::PlantAvailableWater
    "Prognostic crop carbon pool + organ partitioning"
    carbon::Carbon
    "Prognostic crop nitrogen pool + organ partitioning"
    nitrogen::Nitrogen
end

# The default is a well-watered crop (β=1). To make photosynthesis respond to soil water, pass
# `plant_available_water = FieldCapacityLimitedPAW(NF)` AND configure the soil with a clay-bearing
# texture — the default pure-sand texture gives field_capacity == wilting_point, so β is undefined.
function CropVegetation(
        ::Type{NF};
        phenology_dynamics = CropPhenologyDynamics(NF),
        phenology = CropPhenology(NF),
        stomatal_conductance = CropStomatalConductance(NF),
        photosynthesis = CropPhotosynthesis(NF),
        root_distribution = CropRootDistribution(NF),
        plant_available_water = nothing,
        carbon = CropCarbon(NF),
        nitrogen = CropNitrogen(NF),
    ) where {NF}
    return CropVegetation(
        phenology_dynamics, phenology, stomatal_conductance, photosynthesis,
        root_distribution, plant_available_water, carbon, nitrogen,
    )
end

# Skip the plant-available-water computation when it is not configured (β keeps its default of 1).
@inline compute_plant_available_water!(state, grid, ::Nothing, soil) = nothing
@inline compute_plant_available_water!(state, grid, paw, soil) = compute_auxiliary!(state, grid, paw, soil)

"""
    $(TYPEDSIGNATURES)

Compute the crop vegetation auxiliaries in dependency order: soil-moisture stress β (if configured) →
heat-unit fraction → LAI → stomatal conductance (λ) → photosynthesis (GPP).
"""
function Terrarium.compute_auxiliary!(
        state, grid, veg::CropVegetation,
        constants::Terrarium.PhysicalConstants, atmos::Terrarium.AbstractAtmosphere,
        soil::Terrarium.AbstractSoil, args...,
    )
    compute_plant_available_water!(state, grid, veg.plant_available_water, soil)
    compute_auxiliary!(state, grid, veg.phenology_dynamics)
    compute_auxiliary!(state, grid, veg.phenology)
    compute_auxiliary!(state, grid, veg.stomatal_conductance, veg.photosynthesis, constants, atmos)
    compute_auxiliary!(state, grid, veg.photosynthesis, veg.stomatal_conductance, constants, atmos)
    compute_auxiliary!(state, grid, veg.carbon)     # organ carbon partitioning + NPP (needs GPP)
    compute_auxiliary!(state, grid, veg.nitrogen)   # organ nitrogen partitioning (needs organ carbon)
    return nothing
end

"""
    $(TYPEDSIGNATURES)

Initialize the crop vegetation. The `nitrogen_limitation` factor is a lagged auxiliary consumed by
photosynthesis before the nitrogen pool has run, so it is seeded to 1 (no limitation) to avoid
zeroing the first assimilation step.
"""
function Terrarium.initialize!(state, grid, veg::CropVegetation, args...)
    set!(state.nitrogen_limitation, one(eltype(state.nitrogen_limitation)))
    return nothing
end

""" $(TYPEDSIGNATURES) Integrate the prognostic heat units, crop biomass, and crop nitrogen. """
function Terrarium.compute_tendencies!(state, grid, veg::CropVegetation, args...)
    compute_tendencies!(state, grid, veg.phenology_dynamics)
    compute_tendencies!(state, grid, veg.carbon)
    compute_tendencies!(state, grid, veg.nitrogen)
    return nothing
end
