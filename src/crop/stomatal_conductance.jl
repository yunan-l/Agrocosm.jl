# Crop stomatal conductance: supplies the leaf-to-air CO₂ ratio (λ) and canopy water
# conductance that couple crop photosynthesis to the surface water balance.
#
# This continuous first version follows LPJmL's optimal internal:ambient CO₂ ratio: a
# well-watered crop operates near λ = LAMBDA_OPT (0.8), and water limitation reduces λ toward
# a stressed floor in proportion to the soil-moisture limiting factor β. The full LPJmL λ
# water-coupling solver (finding λ where diffusive CO₂ supply equals the biochemical demand via
# the 30-step `lpj_bisect`) is a documented refinement — see
# docs/dev/2026-07/2026-07-23_PHASE3_photosynthesis_design.md. Unlike Terrarium's
# `MedlynStomatalConductance`, this dispatches on any `AbstractPhotosynthesis`, so it pairs with
# `CropPhotosynthesis`.

"""
    $(TYPEDEF)

Crop stomatal conductance. Produces the leaf-to-air CO₂ ratio `λc` and the canopy water
conductance from the soil-moisture limiting factor β and the leaf area index.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropStomatalConductance{NF} <: Terrarium.AbstractStomatalConductance{NF}
    "Optimal (well-watered) internal:ambient CO₂ ratio (LPJmL LAMBDA_OPT)"
    @param λ_opt::NF = 0.8 (bounds = UnitInterval,)
    "Water-stressed lower bound on the internal:ambient CO₂ ratio"
    @param λ_min::NF = 0.2 (bounds = UnitInterval,)
    "Minimum canopy water conductance"
    @param g_min::NF = 1.0e-4 (units = u"m/s", bounds = Positive)
    "Maximum canopy water conductance (well-watered, full canopy)"
    @param g_max::NF = 1.0e-2 (units = u"m/s", bounds = Positive)
    "Canopy light-extinction coefficient"
    @param k_ext::NF = 0.5 (bounds = Positive,)
end

CropStomatalConductance(::Type{NF}; kwargs...) where {NF} = CropStomatalConductance{NF}(; kwargs...)

Terrarium.variables(::CropStomatalConductance) = (
    Terrarium.auxiliary(:canopy_water_conductance, XY(), units = u"m/s"),
    Terrarium.auxiliary(:leaf_to_air_co2_ratio, XY()),
)

@inline Base.@propagate_inbounds Terrarium.stomatal_conductance(i, j, grid, fields, ::CropStomatalConductance) =
    fields.canopy_water_conductance[i, j]

# ---- scalar primitives (Level III) --------------------------------------------------------

"""$(TYPEDSIGNATURES) Leaf-to-air CO₂ ratio: LAMBDA_OPT reduced toward `λ_min` under water stress."""
@inline function compute_leaf_to_air_co2_ratio(stomcond::CropStomatalConductance{NF}, β::NF) where {NF}
    return stomcond.λ_min + (stomcond.λ_opt - stomcond.λ_min) * β
end

"""$(TYPEDSIGNATURES) Canopy water conductance scaled by canopy cover (Lambert–Beer) and β."""
@inline function compute_canopy_conductance(stomcond::CropStomatalConductance{NF}, LAI::NF, β::NF) where {NF}
    canopy_cover = NF(1) - exp(-stomcond.k_ext * LAI)
    return stomcond.g_min + (stomcond.g_max - stomcond.g_min) * canopy_cover * β
end

# ---- kernel functions (Level II) ----------------------------------------------------------

"""$(TYPEDSIGNATURES) Canopy water conductance and leaf-to-air CO₂ ratio at grid point `(i, j)`."""
Base.@propagate_inbounds function compute_stomatal_conductance(
        i, j, grid, fields, stomcond::CropStomatalConductance,
        photo::Terrarium.AbstractPhotosynthesis, constants::Terrarium.PhysicalConstants,
        atmos::Terrarium.AbstractAtmosphere, args...,
    )
    LAI = fields.leaf_area_index[i, j]
    # β is a limiting factor by definition bounded to [0, 1]; clamp defensively so an
    # upstream out-of-range value cannot drive the conductance or λ outside their ranges.
    β = clamp(fields.soil_moisture_limiting_factor[i, j], zero(LAI), one(LAI))
    g_stm = compute_canopy_conductance(stomcond, LAI, β)
    λc = compute_leaf_to_air_co2_ratio(stomcond, β)
    return g_stm, λc
end

"""$(TYPEDSIGNATURES) Store [`compute_stomatal_conductance`](@ref) outputs in `out`."""
Base.@propagate_inbounds function compute_stomatal_conductance!(
        out, i, j, grid, fields, stomcond::CropStomatalConductance,
        photo::Terrarium.AbstractPhotosynthesis, constants::Terrarium.PhysicalConstants,
        atmos::Terrarium.AbstractAtmosphere, args...,
    )
    g_stm, λc = compute_stomatal_conductance(i, j, grid, fields, stomcond, photo, constants, atmos, args...)
    out.canopy_water_conductance[i, j, 1] = g_stm
    out.leaf_to_air_co2_ratio[i, j, 1] = λc
    return out
end

# ---- interface methods (Level I) ----------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(
        state, grid, stomcond::CropStomatalConductance,
        photo::Terrarium.AbstractPhotosynthesis, constants::Terrarium.PhysicalConstants,
        atmos::Terrarium.AbstractAtmosphere, args...,
    )
    out = Terrarium.auxiliary_fields(state, stomcond)
    # `photo` declares leaf_area_index and soil_moisture_limiting_factor as inputs.
    fields = get_fields(state, stomcond, photo, constants, atmos; except = out)
    launch!(grid, XY, compute_stomatal_conductance_kernel!, out, fields, stomcond, photo, constants, atmos)
    return nothing
end

@kernel inbounds = true function compute_stomatal_conductance_kernel!(out, grid, fields, stomcond::CropStomatalConductance, photo, constants, atmos)
    i, j = @index(Global, NTuple)
    compute_stomatal_conductance!(out, i, j, grid, fields, stomcond, photo, constants, atmos)
end
