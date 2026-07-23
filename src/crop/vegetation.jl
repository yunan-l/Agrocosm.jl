# A minimal crop vegetation model wiring the crop processes into Terrarium's `AbstractVegetation`
# interface, so it slots into a `LandModel` alongside Terrarium's soil, surface, and atmosphere. The
# crop leaf area index is driven by the prognostic phenological heat units (not a carbon-pool
# equilibrium): heat units accumulate with temperature, set the heat-unit fraction, which drives the
# LAI trajectory, which feeds photosynthesis. This is the Phase 5 assembly of the ported crop
# primitives; the multi-organ carbon/nitrogen pools and soil C–N coupling are added incrementally.

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
    } <: Terrarium.AbstractVegetation{NF}
    "Phenological heat-unit accumulation (prognostic)"
    phenology_dynamics::PhenologyDynamics
    "Heat-unit-driven leaf-area-index trajectory"
    phenology::Phenology
    "Crop stomatal conductance (λ, canopy conductance)"
    stomatal_conductance::StomatalConductance
    "C3/C4 crop photosynthesis"
    photosynthesis::Photosynthesis
end

function CropVegetation(
        ::Type{NF};
        phenology_dynamics = CropPhenologyDynamics(NF),
        phenology = CropPhenology(NF),
        stomatal_conductance = CropStomatalConductance(NF),
        photosynthesis = CropPhotosynthesis(NF),
    ) where {NF}
    return CropVegetation(phenology_dynamics, phenology, stomatal_conductance, photosynthesis)
end

"""
    $(TYPEDSIGNATURES)

Compute the crop vegetation auxiliaries in dependency order: heat-unit fraction → LAI → stomatal
conductance (λ) → photosynthesis (GPP). `soil` and any further arguments are accepted for interface
compatibility but not used by this minimal model.
"""
function Terrarium.compute_auxiliary!(
        state, grid, veg::CropVegetation,
        constants::Terrarium.PhysicalConstants, atmos::Terrarium.AbstractAtmosphere, args...,
    )
    compute_auxiliary!(state, grid, veg.phenology_dynamics)
    compute_auxiliary!(state, grid, veg.phenology)
    compute_auxiliary!(state, grid, veg.stomatal_conductance, veg.photosynthesis, constants, atmos)
    compute_auxiliary!(state, grid, veg.photosynthesis, veg.stomatal_conductance, constants, atmos)
    return nothing
end

""" $(TYPEDSIGNATURES) Integrate the prognostic phenological heat units. """
function Terrarium.compute_tendencies!(state, grid, veg::CropVegetation, args...)
    compute_tendencies!(state, grid, veg.phenology_dynamics)
    return nothing
end
