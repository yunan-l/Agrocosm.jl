# Crop C3/C4 photosynthesis as a continuous-time Terrarium `AbstractPhotosynthesis`.
#
# See docs/dev/2026-07/2026-07-23_PHASE3_photosynthesis_design.md. The C3 pathway is the
# BIOME3/Haxeltine & Prentice (1996) mechanistic scheme shared with Terrarium's
# `LUEPhotosynthesis` (validated by an equivalence test); the C4 pathway drops the
# Γ*/Kc/Ko chain (c₂ ≡ 1) and scales the light term by φ = min(1, λ/λ_mc4). Instantaneous
# rates (gC/m²/s), integrated by the timestepper — the legacy daily/daylength scaling is
# dropped. λ (the intercellular:ambient CO₂ ratio) is supplied as the input
# `leaf_to_air_co2_ratio` by a stomatal-conductance process.

"""Photosynthetic pathway discriminant for [`CropPhotosynthesis`](@ref)."""
abstract type AbstractCropPathway end

"""C3 photosynthetic pathway (e.g. wheat, rice, soybean)."""
struct C3Pathway <: AbstractCropPathway end

"""C4 photosynthetic pathway (e.g. maize, millet, sugarcane)."""
struct C4Pathway <: AbstractCropPathway end

"""
    $(TYPEDEF)

Crop C3/C4 photosynthesis. The `pathway` field selects the biochemistry; parameters default
to the PALADYN/BIOME3 needleleaf values so a C3 instance reproduces Terrarium's
`LUEPhotosynthesis`. Crop functional-type presets (temperature thresholds, quantum
efficiencies) are wired from the CFT registry in Phase 5.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropPhotosynthesis{NF, Path <: AbstractCropPathway} <: Terrarium.AbstractPhotosynthesis{NF}
    "Photosynthetic pathway (`C3Pathway()` or `C4Pathway()`)"
    pathway::Path = C3Pathway()

    "Rubisco specificity factor at 25 °C (C3)"
    @param τ25::NF = 2600.0 (bounds = Positive,)
    "Michaelis-Menten constant for CO₂ at 25 °C (C3)"
    @param Kc25::NF = 30.0 (units = u"Pa", bounds = Positive)
    "Michaelis-Menten constant for O₂ at 25 °C (C3)"
    @param Ko25::NF = 3.0e4 (units = u"Pa", bounds = Positive)
    "Q10 temperature sensitivity of τ"
    @param q10_τ::NF = 0.57 (bounds = Positive,)
    "Q10 temperature sensitivity of Kc"
    @param q10_Kc::NF = 2.1 (bounds = Positive,)
    "Q10 temperature sensitivity of Ko"
    @param q10_Ko::NF = 1.2 (bounds = Positive,)
    "Leaf albedo in the PAR range"
    @param α_leaf::NF = 0.17 (bounds = UnitInterval,)
    "Fraction of PAR assimilated at ecosystem level"
    @param α_a::NF = 0.5 (bounds = UnitInterval,)
    "Intrinsic quantum efficiency of CO₂ uptake, C3"
    @param α_C3::NF = 0.08 (bounds = UnitInterval,)
    "Intrinsic quantum efficiency of CO₂ uptake, C4"
    @param α_C4::NF = 0.053 (bounds = UnitInterval,)
    "PAR photon-to-mole conversion (mol photons per J)"
    @param cq::NF = 4.6e-6 (units = u"mol/J", bounds = Positive, scale = 1.0e-6)
    "Canopy light-extinction coefficient"
    @param k_ext::NF = 0.5 (bounds = Positive,)
    "Upper temperature threshold for CO₂/O₂ specificity"
    @param T_CO2_high::NF = 42.0 (units = u"°C",)
    "Lower temperature threshold for CO₂/O₂ specificity"
    @param T_CO2_low::NF = -4.0 (units = u"°C",)
    "Upper temperature threshold for light-limited photosynthesis"
    @param T_photos_high::NF = 30.0 (units = u"°C",)
    "Lower temperature threshold for light-limited photosynthesis"
    @param T_photos_low::NF = 15.0 (units = u"°C",)
    "Co-limitation shape parameter"
    @param θ_r::NF = 0.7 (bounds = UnitInterval,)
    "High-temperature photosynthesis cutoff, C3 (LPJmL tmc3)"
    @param T_cutoff_C3::NF = 45.0 (units = u"°C",)
    "High-temperature photosynthesis cutoff, C4 (LPJmL tmc4)"
    @param T_cutoff_C4::NF = 55.0 (units = u"°C",)
    "C4 nominal internal:ambient CO₂ ratio (φ saturation)"
    @param λ_mc4::NF = 0.4 (bounds = UnitInterval,)
end

CropPhotosynthesis(::Type{NF}; pathway::AbstractCropPathway = C3Pathway(), kwargs...) where {NF} =
    CropPhotosynthesis{NF, typeof(pathway)}(; pathway, kwargs...)

Terrarium.variables(::CropPhotosynthesis{NF}) where {NF} = (
    Terrarium.auxiliary(:net_assimilation, XY(), units = u"g/m^2/s"),
    Terrarium.auxiliary(:leaf_respiration, XY(), units = u"g/m^2/s"),
    Terrarium.auxiliary(:gross_primary_production, XY(), units = u"kg/m^2/s"),
    Terrarium.input(:soil_moisture_limiting_factor, XY(), default = NF(1)),
    Terrarium.input(:leaf_area_index, XY()),
)

# ---- scalar primitives (Level III) --------------------------------------------------------

"""$(TYPEDSIGNATURES) Net PAR reaching the canopy (mol/m²/s)."""
@inline function compute_par(photo::CropPhotosynthesis{NF}, swdown::NF) where {NF}
    return NF(0.5) * swdown * (NF(1) - photo.α_leaf) * photo.cq
end

"""$(TYPEDSIGNATURES) Absorbed PAR (mol/m²/s), Lambert–Beer over the canopy."""
@inline function compute_apar(photo::CropPhotosynthesis{NF}, par::NF, LAI::NF) where {NF}
    return photo.α_a * par * (NF(1) - exp(-photo.k_ext * LAI))
end

"""$(TYPEDSIGNATURES) Double-sigmoid temperature stress with a pathway-specific high-T cutoff."""
@inline function compute_temperature_stress(photo::CropPhotosynthesis{NF}, path::AbstractCropPathway, T_air::NF) where {NF}
    k1 = NF(2) * log(NF(1) / NF(0.99) - NF(1)) / (photo.T_CO2_low - photo.T_photos_low)
    k2 = NF(0.5) * (photo.T_CO2_low + photo.T_photos_low)
    k3 = log(NF(0.99) / NF(0.01)) / (photo.T_CO2_high - photo.T_photos_high)
    cutoff = high_temperature_cutoff(photo, path)
    if photo.T_CO2_low < T_air < photo.T_CO2_high && T_air ≤ cutoff
        low = NF(1) / (NF(1) + exp(k1 * (k2 - T_air)))
        high = NF(1) - NF(0.01) * exp(k3 * (T_air - photo.T_photos_high))
        return low * high
    else
        return zero(NF)
    end
end

@inline high_temperature_cutoff(photo::CropPhotosynthesis, ::C3Pathway) = photo.T_cutoff_C3
@inline high_temperature_cutoff(photo::CropPhotosynthesis, ::C4Pathway) = photo.T_cutoff_C4

"""
    $(TYPEDSIGNATURES)

Pathway-specific light-limited factor `c₁` (gC/mol), Rubisco-limited factor `c₂`, and the
maximum carboxylation rate `Vc_max` (gC/m²/s), from the intercellular CO₂ pressure `pres_i`,
temperature stress, PAR, and the Q10 kinetics (C3 only).
"""
# Note: following BIOME3/LUE, the coordination hypothesis sets Vc_max from *absorbed* PAR
# (APAR), not incident PAR — so JC = c₂·Vc_max coordinates with the light-limited JE = c₁·APAR.
@inline function compute_assimilation_terms(
        photo::CropPhotosynthesis{NF}, ::C3Pathway, cmass::NF, T_air::NF, T_stress::NF,
        apar::NF, pres_i::NF, pres_O2::NF, λc::NF,
    ) where {NF}
    τ = photo.τ25 * photo.q10_τ^((T_air - NF(25)) * NF(0.1))
    Kc = photo.Kc25 * photo.q10_Kc^((T_air - NF(25)) * NF(0.1))
    Ko = photo.Ko25 * photo.q10_Ko^((T_air - NF(25)) * NF(0.1))
    Γ_star = pres_O2 / (NF(2) * τ)
    c_1 = photo.α_C3 * T_stress * cmass * (pres_i - Γ_star) / (pres_i + NF(2) * Γ_star)
    c_2 = (pres_i - Γ_star) / (pres_i + Kc * (NF(1) + pres_O2 / Ko))
    Vc_max = c_1 * apar * (pres_i + Kc * (NF(1) + pres_O2 / Ko)) / (pres_i - Γ_star)
    return c_1, c_2, Vc_max
end

@inline function compute_assimilation_terms(
        photo::CropPhotosynthesis{NF}, ::C4Pathway, cmass::NF, T_air::NF, T_stress::NF,
        apar::NF, pres_i::NF, pres_O2::NF, λc::NF,
    ) where {NF}
    # C4: the CO₂-concentrating mechanism removes the Rubisco/oxygenation limitation (c₂ ≡ 1);
    # the light term saturates through φ = min(1, λ/λ_mc4) in the intercellular:ambient ratio λ.
    φ = min(NF(1), λc / photo.λ_mc4)
    c_1 = photo.α_C4 * T_stress * cmass * φ
    c_2 = one(NF)
    Vc_max = c_1 * apar
    return c_1, c_2, Vc_max
end

"""$(TYPEDSIGNATURES) θ-form co-limitation of light (`JE`) and Rubisco (`JC`) rates."""
@inline function co_limit(photo::CropPhotosynthesis{NF}, JE::NF, JC::NF) where {NF}
    discriminant = max(zero(NF), (JE + JC)^2 - NF(4) * photo.θ_r * JE * JC)
    return (JE + JC - sqrt(discriminant)) / (NF(2) * photo.θ_r)
end

"""$(TYPEDSIGNATURES) α coefficient for leaf/maintenance respiration `Rd = α·Vc_max·β`."""
@inline respiration_coefficient(photo::CropPhotosynthesis, ::C3Pathway) = photo.α_C3
@inline respiration_coefficient(photo::CropPhotosynthesis, ::C4Pathway) = photo.α_C4

"""
    $(TYPEDSIGNATURES)

Leaf respiration `Rd` and net assimilation `An` (gC/m²/s) at a single point, for the given
air temperature, incoming shortwave, pressure, CO₂ (ppm), LAI, intercellular:ambient CO₂
ratio `λc`, and soil-moisture limiting factor `β`.
"""
function compute_respiration_assimilation(
        photo::CropPhotosynthesis{NF}, cmass::NF,
        T_air::NF, swdown::NF, pres::NF, co2::NF, LAI::NF, λc::NF, β::NF,
    ) where {NF}
    path = photo.pathway
    pres_O2 = Terrarium.partial_pressure_O2(pres)
    pres_a = Terrarium.partial_pressure_CO2(pres, co2)
    if swdown > zero(NF) && T_air > NF(-3) && LAI > zero(NF)
        par = compute_par(photo, swdown)
        apar = compute_apar(photo, par, LAI)
        pres_i = λc * pres_a
        T_stress = compute_temperature_stress(photo, path, T_air)
        c_1, c_2, Vc_max = compute_assimilation_terms(photo, path, cmass, T_air, T_stress, apar, pres_i, pres_O2, λc)
        Rd = respiration_coefficient(photo, path) * Vc_max * β
        JE = c_1 * apar
        JC = c_2 * Vc_max
        Ag = co_limit(photo, JE, JC) * β
        An = Ag - Rd
        return Rd, An
    else
        return zero(NF), zero(NF)
    end
end

"""$(TYPEDSIGNATURES) Gross primary production (kgC/m²/s) from net assimilation (gC/m²/s)."""
@inline compute_gpp(::CropPhotosynthesis{NF}, An::NF) where {NF} = An * NF(1.0e-3)

# ---- kernel functions (Level II) ----------------------------------------------------------

"""$(TYPEDSIGNATURES) Leaf respiration, net assimilation, and GPP at grid point `(i, j)`."""
Base.@propagate_inbounds function compute_photosynthesis(
        i, j, grid, fields, photo::CropPhotosynthesis,
        constants::Terrarium.PhysicalConstants, atmos::Terrarium.AbstractAtmosphere,
    )
    T_air = Terrarium.air_temperature(i, j, grid, fields, atmos)
    pres = Terrarium.air_pressure(i, j, grid, fields, atmos)
    swdown = Terrarium.shortwave_down(i, j, grid, fields, atmos)
    co2 = fields.CO2[i, j]
    β = fields.soil_moisture_limiting_factor[i, j]
    LAI = fields.leaf_area_index[i, j]
    λc = fields.leaf_to_air_co2_ratio[i, j]
    cmass = constants.material.atomic_weight_carbon
    Rd, An = compute_respiration_assimilation(photo, cmass, T_air, swdown, pres, co2, LAI, λc, β)
    GPP = compute_gpp(photo, An)
    return Rd, An, GPP
end

"""$(TYPEDSIGNATURES) Store [`compute_photosynthesis`](@ref) outputs in `out`."""
Base.@propagate_inbounds function compute_photosynthesis!(
        out, i, j, grid, fields, photo::CropPhotosynthesis,
        constants::Terrarium.PhysicalConstants, atmos::Terrarium.AbstractAtmosphere,
    )
    Rd, An, GPP = compute_photosynthesis(i, j, grid, fields, photo, constants, atmos)
    out.leaf_respiration[i, j, 1] = Rd
    out.net_assimilation[i, j, 1] = An
    out.gross_primary_production[i, j, 1] = GPP
    return out
end

# ---- interface methods (Level I) ----------------------------------------------------------

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(
        state, grid, photo::CropPhotosynthesis,
        stomcond::Terrarium.AbstractStomatalConductance,
        constants::Terrarium.PhysicalConstants,
        atmos::Terrarium.AbstractAtmosphere, args...,
    )
    out = Terrarium.auxiliary_fields(state, photo)
    fields = get_fields(state, photo, stomcond, atmos; except = out)
    launch!(grid, XY, compute_photosynthesis_kernel!, out, fields, photo, constants, atmos)
    return nothing
end

@kernel inbounds = true function compute_photosynthesis_kernel!(out, grid, fields, photo::CropPhotosynthesis, constants, atmos)
    i, j = @index(Global, NTuple)
    compute_photosynthesis!(out, i, j, grid, fields, photo, constants, atmos)
end
