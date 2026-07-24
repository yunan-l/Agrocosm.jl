# Crop soil carbon biogeochemistry — a Terrarium `AbstractSoilBiogeochemistry` with prognostic soil
# carbon pools (litter, fast, slow; kgC/m³) that decompose by first-order rates modulated by the
# environmental decomposition response (soil temperature × moisture). Decomposed litter is routed to
# the fast/slow pools and the atmosphere; decomposed fast/slow carbon is respired. The live fast+slow
# density feeds `density_soc`, so the soil organic fraction (and hence porosity and thermal/hydraulic
# properties) responds to the carbon dynamics. Rates are per-day (LPJmL) applied per second.
#
# This replaces the constant `ConstantSoilCarbonDensity` in the soil `biogeochem` slot with dynamic
# pools. Nitrogen transforms (nitrification/denitrification/mineralization — already ported as tested
# primitives) and the crop-litterfall input into the litter pool are the next coupling steps.

"""
    $(TYPEDEF)

Prognostic crop soil carbon biogeochemistry (litter/fast/slow decomposition).

Properties:
$(TYPEDFIELDS)
"""
@kwdef struct CropSoilBiogeochemistry{NF} <: Terrarium.AbstractSoilBiogeochemistry{NF}
    "Carbon pool decomposition rates + litter routing fractions"
    carbon::CropSoilCarbon{NF} = CropSoilCarbon(NF)
    "Environmental (temperature × moisture) decomposition response"
    response::CropSoilDecompositionResponse{NF} = CropSoilDecompositionResponse(NF)
    "Nitrification (NH₄ → NO₃) parameters"
    nitrification::CropNitrification{NF} = CropNitrification(NF)
    "Denitrification (NO₃ → gas) parameters"
    denitrification::CropDenitrification{NF} = CropDenitrification(NF)
    "Litter decomposition rate at 10 °C"
    k_litter::NF = 0.5 / 365
    "Soil organic-matter C:N ratio governing net mineralization"
    soil_cn_ratio::NF = 15.0
    "Soil pH (for the nitrification response)"
    soil_ph::NF = 6.5
    "Pure organic-matter density (for the organic solid fraction)"
    ρ_org::NF = 1300.0
    "Initial fast-pool carbon density"
    initial_fast_carbon::NF = 5.0
    "Initial slow-pool carbon density"
    initial_slow_carbon::NF = 20.0
    "Initial litter carbon density"
    initial_litter_carbon::NF = 1.0
    "Initial soil ammonium density"
    initial_ammonium::NF = 0.05
    "Initial soil nitrate density"
    initial_nitrate::NF = 0.05
end

CropSoilBiogeochemistry(::Type{NF}; kwargs...) where {NF} = CropSoilBiogeochemistry{NF}(; kwargs...)

Terrarium.variables(::CropSoilBiogeochemistry{NF}) where {NF} = (
    Terrarium.prognostic(:litter_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:fast_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:slow_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:soil_ammonium, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:soil_nitrate, XYZ(), units = u"kg/m^3"),
    Terrarium.auxiliary(:decomposition_response, XYZ()),
    Terrarium.auxiliary(:heterotrophic_respiration, XYZ(), units = u"kg/m^3/s"),
    Terrarium.auxiliary(:net_mineralization, XYZ(), units = u"kg/m^3/s"),
    Terrarium.input(:temperature, XYZ(), default = NF(5), units = u"°C"),
    Terrarium.input(:saturation_water_ice, XYZ(), default = NF(0.5)),
    # Crop coupling: 0D per-area fluxes distributed over the root zone, and the root fraction.
    Terrarium.input(:root_fraction, XYZ(), default = zero(NF)),
    Terrarium.input(:crop_litterfall_carbon, XY(), default = zero(NF), units = u"kg/m^2/s"),
    Terrarium.input(:crop_litterfall_nitrogen, XY(), default = zero(NF), units = u"kg/m^2/s"),
    Terrarium.input(:crop_nitrogen_uptake, XY(), default = zero(NF), units = u"kg/m^2/s"),
)

Terrarium.density_pure_soc(bgc::CropSoilBiogeochemistry) = bgc.ρ_org

"""$(TYPEDSIGNATURES) Organic carbon density (kgC/m³) as the live fast + slow soil carbon."""
Base.@propagate_inbounds Terrarium.density_soc(i, j, k, grid, fields, ::CropSoilBiogeochemistry) =
    fields.fast_carbon[i, j, k] + fields.slow_carbon[i, j, k]

"""
    $(TYPEDSIGNATURES)

Per-second pool tendencies `(d_litter, d_fast, d_slow)` and heterotrophic respiration, from the
current pools and the environmental decomposition `response`. First-order decay `λ_x = k_x·response`
(per day) is applied per second; decomposed litter is routed to fast/slow/atmosphere.
"""
@inline function soil_carbon_tendencies(bgc::CropSoilBiogeochemistry{NF}, litter::NF, fast::NF, slow::NF, response::NF) where {NF}
    per_second = response / Terrarium.seconds_per_day(NF)
    decomposed_litter = bgc.k_litter * per_second * max(zero(NF), litter)
    decomposed_fast = bgc.carbon.k_fast * per_second * max(zero(NF), fast)
    decomposed_slow = bgc.carbon.k_slow * per_second * max(zero(NF), slow)
    to_fast, to_slow, to_atmosphere = route_litter_carbon(bgc.carbon, decomposed_litter)
    d_litter = -decomposed_litter
    d_fast = to_fast - decomposed_fast
    d_slow = to_slow - decomposed_slow
    heterotrophic_respiration = to_atmosphere + decomposed_fast + decomposed_slow
    return d_litter, d_fast, d_slow, heterotrophic_respiration
end

"""
    $(TYPEDSIGNATURES)

Per-second mineral-nitrogen tendencies `(d_ammonium, d_nitrate)` from the current NH₄/NO₃ pools, the
net mineralization (N released by the respired carbon, kgN/m³/s), the soil temperature, the
water-filled pore space, and the organic (fast + slow) carbon. Mineralization feeds NH₄; nitrification
moves NH₄ → NO₃ (minus the N₂O loss); denitrification removes NO₃.
"""
@inline function soil_nitrogen_tendencies(
        bgc::CropSoilBiogeochemistry{NF}, ammonium::NF, nitrate::NF, mineralization::NF,
        temperature::NF, water_filled_pore_space::NF, organic_carbon::NF,
    ) where {NF}
    per_second = one(NF) / Terrarium.seconds_per_day(NF)
    gross_nit, n2o_nit = gross_nitrification(bgc.nitrification, max(zero(NF), ammonium), water_filled_pore_space, temperature, bgc.soil_ph)
    gross_denit, _n2o, _n2 = gross_denitrification(bgc.denitrification, max(zero(NF), nitrate), temperature, water_filled_pore_space, organic_carbon)
    nitrification_rate = gross_nit * per_second
    n2o_nitrification_rate = n2o_nit * per_second
    denitrification_rate = gross_denit * per_second
    d_ammonium = mineralization - nitrification_rate
    d_nitrate = (nitrification_rate - n2o_nitrification_rate) - denitrification_rate
    return d_ammonium, d_nitrate
end

# ---- interface methods --------------------------------------------------------------------

""" $(TYPEDSIGNATURES) Seed the soil carbon and mineral-nitrogen pools. """
function Terrarium.initialize!(state, grid, bgc::CropSoilBiogeochemistry, args...)
    set!(state.litter_carbon, bgc.initial_litter_carbon)
    set!(state.fast_carbon, bgc.initial_fast_carbon)
    set!(state.slow_carbon, bgc.initial_slow_carbon)
    set!(state.soil_ammonium, bgc.initial_ammonium)
    set!(state.soil_nitrate, bgc.initial_nitrate)
    return nothing
end

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_auxiliary!(state, grid, bgc::CropSoilBiogeochemistry, args...)
    out = Terrarium.auxiliary_fields(state, bgc)
    fields = get_fields(state, bgc; except = out)
    launch!(grid, XYZ, compute_soil_bgc_auxiliary_kernel!, out, fields, bgc)
    return nothing
end

""" $(TYPEDSIGNATURES) """
function Terrarium.compute_tendencies!(state, grid, bgc::CropSoilBiogeochemistry, args...)
    out = Terrarium.tendency_fields(state, bgc)
    fields = get_fields(state, bgc)
    launch!(grid, XYZ, compute_soil_bgc_tendency_kernel!, out, fields, bgc)
    return nothing
end

@kernel inbounds = true function compute_soil_bgc_auxiliary_kernel!(out, grid, fields, bgc::CropSoilBiogeochemistry)
    i, j, k = @index(Global, NTuple)
    resp = soil_decomposition_response(bgc.response, fields.temperature[i, j, k], fields.saturation_water_ice[i, j, k])
    out.decomposition_response[i, j, k] = resp
    _dl, _df, _ds, het = soil_carbon_tendencies(
        bgc, fields.litter_carbon[i, j, k], fields.fast_carbon[i, j, k], fields.slow_carbon[i, j, k], resp,
    )
    out.heterotrophic_respiration[i, j, k] = het
    # Net mineralization: nitrogen released by the respired carbon at the soil C:N ratio.
    out.net_mineralization[i, j, k] = het / bgc.soil_cn_ratio
end

@kernel inbounds = true function compute_soil_bgc_tendency_kernel!(out, grid, fields, bgc::CropSoilBiogeochemistry)
    i, j, k = @index(Global, NTuple)
    NF = eltype(out.litter_carbon)
    resp = fields.decomposition_response[i, j, k]
    d_litter, d_fast, d_slow, _het = soil_carbon_tendencies(
        bgc, fields.litter_carbon[i, j, k], fields.fast_carbon[i, j, k], fields.slow_carbon[i, j, k], resp,
    )
    ammonium = fields.soil_ammonium[i, j, k]
    nitrate = fields.soil_nitrate[i, j, k]
    organic_carbon = fields.fast_carbon[i, j, k] + fields.slow_carbon[i, j, k]
    d_ammonium, d_nitrate = soil_nitrogen_tendencies(
        bgc, ammonium, nitrate, fields.net_mineralization[i, j, k],
        fields.temperature[i, j, k], fields.saturation_water_ice[i, j, k], organic_carbon,
    )

    # Crop coupling: distribute the 0D per-area crop fluxes over the root zone as per-volume rates
    # (÷ layer thickness); the root fraction sums to unity over the column, so mass is conserved.
    field_grid = get_field_grid(grid)
    per_volume = fields.root_fraction[i, j, k] / Δzᵃᵃᶜ(i, j, k, field_grid)
    litterfall_carbon = fields.crop_litterfall_carbon[i, j] * per_volume
    litterfall_nitrogen = fields.crop_litterfall_nitrogen[i, j] * per_volume
    uptake = fields.crop_nitrogen_uptake[i, j] * per_volume
    # Split the crop uptake between the ammonium and nitrate pools by their share.
    total_mineral = max(ammonium + nitrate, eps(NF))
    uptake_ammonium = uptake * ammonium / total_mineral
    uptake_nitrate = uptake * nitrate / total_mineral

    out.litter_carbon[i, j, k] = d_litter + litterfall_carbon
    out.fast_carbon[i, j, k] = d_fast
    out.slow_carbon[i, j, k] = d_slow
    # Litterfall nitrogen mineralizes into ammonium; crop uptake draws down both mineral pools.
    out.soil_ammonium[i, j, k] = d_ammonium + litterfall_nitrogen - uptake_ammonium
    out.soil_nitrate[i, j, k] = d_nitrate - uptake_nitrate
end
