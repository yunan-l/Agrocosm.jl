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
    "Litter mineralization/immobilization (microbial N demand of the decomposing litter)"
    mineralization::CropNitrogenMineralization{NF} = CropNitrogenMineralization(NF)
    "Ammonia volatilization (NH₃ loss from the top soil layer)"
    volatilization::CropVolatilization{NF} = CropVolatilization(NF)
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
    "Initial litter nitrogen density (litter at the soil C:N: initial_litter_carbon / soil_cn_ratio)"
    initial_litter_nitrogen::NF = 1.0 / 15.0
    "Initial soil ammonium density"
    initial_ammonium::NF = 0.05
    "Initial soil nitrate density"
    initial_nitrate::NF = 0.05
end

CropSoilBiogeochemistry(::Type{NF}; kwargs...) where {NF} = CropSoilBiogeochemistry{NF}(; kwargs...)

Terrarium.variables(::CropSoilBiogeochemistry{NF}) where {NF} = (
    Terrarium.prognostic(:litter_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:litter_nitrogen, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:fast_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:slow_carbon, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:soil_ammonium, XYZ(), units = u"kg/m^3"),
    Terrarium.prognostic(:soil_nitrate, XYZ(), units = u"kg/m^3"),
    Terrarium.auxiliary(:decomposition_response, XYZ()),
    Terrarium.auxiliary(:heterotrophic_respiration, XYZ(), units = u"kg/m^3/s"),
    Terrarium.auxiliary(:net_mineralization, XYZ(), units = u"kg/m^3/s"),
    Terrarium.auxiliary(:litter_nitrogen_release, XYZ(), units = u"kg/m^3/s"),
    Terrarium.auxiliary(:ammonia_volatilization, XYZ(), units = u"kg/m^3/s"),
    Terrarium.input(:temperature, XYZ(), default = NF(5), units = u"°C"),
    Terrarium.input(:saturation_water_ice, XYZ(), default = NF(0.5)),
    # Crop coupling: 0D per-area fluxes distributed over the root zone, and the root fraction.
    Terrarium.input(:root_fraction, XYZ(), default = zero(NF)),
    Terrarium.input(:crop_litterfall_carbon, XY(), default = zero(NF), units = u"kg/m^2/s"),
    Terrarium.input(:crop_litterfall_nitrogen, XY(), default = zero(NF), units = u"kg/m^2/s"),
    Terrarium.input(:crop_nitrogen_uptake, XY(), default = zero(NF), units = u"kg/m^2/s"),
    # Volatilization forcing: near-surface air temperature and wind speed drive the NH₃ mass transfer.
    Terrarium.input(:air_temperature, XY(), default = NF(10), units = u"°C"),
    Terrarium.input(:windspeed, XY(), default = NF(0.1), units = u"m/s"),
    # Fertilizer: continuous 0D per-area mineral-N application fluxes (see `CropFertilization`).
    Terrarium.input(:fertilizer_ammonium_flux, XY(), default = zero(NF), units = u"kg/m^2/s"),
    Terrarium.input(:fertilizer_nitrate_flux, XY(), default = zero(NF), units = u"kg/m^2/s"),
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
    decomposed_litter, decomposed_fast, decomposed_slow, to_fast, to_slow, to_atmosphere =
        soil_carbon_decomposition(bgc, litter, fast, slow, response)
    d_litter = -decomposed_litter
    d_fast = to_fast - decomposed_fast
    d_slow = to_slow - decomposed_slow
    heterotrophic_respiration = to_atmosphere + decomposed_fast + decomposed_slow
    return d_litter, d_fast, d_slow, heterotrophic_respiration
end

"""
    $(TYPEDSIGNATURES)

Per-second carbon-decomposition fluxes `(decomposed_litter, decomposed_fast, decomposed_slow, to_fast,
to_slow, to_atmosphere)` (kgC/m³/s) from the current pools and the environmental `response`. First-order
decay `λ_x = k_x·response` (per day) applied per second; decomposed litter routes to fast/slow/atmosphere.
Shared by the carbon tendencies and the litter-nitrogen accounting.
"""
@inline function soil_carbon_decomposition(bgc::CropSoilBiogeochemistry{NF}, litter::NF, fast::NF, slow::NF, response::NF) where {NF}
    per_second = response / Terrarium.seconds_per_day(NF)
    decomposed_litter = bgc.k_litter * per_second * max(zero(NF), litter)
    decomposed_fast = bgc.carbon.k_fast * per_second * max(zero(NF), fast)
    decomposed_slow = bgc.carbon.k_slow * per_second * max(zero(NF), slow)
    to_fast, to_slow, to_atmosphere = route_litter_carbon(bgc.carbon, decomposed_litter)
    return decomposed_litter, decomposed_fast, decomposed_slow, to_fast, to_slow, to_atmosphere
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

"""
    $(TYPEDSIGNATURES)

Per-volume nitrogen exchange (kgN/m³/s) between the decomposing litter and the mineral pool. As litter
carbon decomposes it liberates nitrogen at the litter's *own* C:N (`n_release`); the carbon routed into
the fast/slow soil organic matter is built at the soil C:N and consumes `(to_fast + to_slow)/soil_cn`
of nitrogen. The surplus `n_release − demand` is mineralized to ammonium when positive; when negative
the deficit is *immobilized* from the mineral pool, Michaelis–Menten-limited by the available ammonium
concentration. Returns `(n_release, litter_to_mineral)`, with `litter_to_mineral` the net litter→ammonium
flux (negative = immobilization drawn from the pool).
"""
@inline function litter_nitrogen_exchange(
        bgc::CropSoilBiogeochemistry{NF}, litter_carbon::NF, litter_nitrogen::NF,
        decomposed_litter::NF, to_fast::NF, to_slow::NF, available_ammonium::NF, layer_depth::NF,
    ) where {NF}
    n_release = decomposed_litter * max(zero(NF), litter_nitrogen) / max(litter_carbon, eps(NF))
    demand = (to_fast + to_slow) / bgc.soil_cn_ratio
    surplus = n_release - demand
    # Immobilization when the litter is nitrogen-poor: the deficit rate is limited by the available
    # mineral-N concentration (per-area gN/m² and layer depth in mm, matching the LPJmL primitive).
    available_perarea = max(zero(NF), available_ammonium) * layer_depth * NF(1000)
    limitation = immobilization_limitation(bgc.mineralization, available_perarea, layer_depth * NF(1000))
    immobilized = max(zero(NF), -surplus) * limitation
    litter_to_mineral = ifelse(surplus >= zero(NF), surplus, -immobilized)
    return n_release, litter_to_mineral
end

# ---- interface methods --------------------------------------------------------------------

""" $(TYPEDSIGNATURES) Seed the soil carbon and mineral-nitrogen pools. """
function Terrarium.initialize!(state, grid, bgc::CropSoilBiogeochemistry, args...)
    set!(state.litter_carbon, bgc.initial_litter_carbon)
    set!(state.litter_nitrogen, bgc.initial_litter_nitrogen)
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
    NF = eltype(out.net_mineralization)
    field_grid = get_field_grid(grid)
    resp = soil_decomposition_response(bgc.response, fields.temperature[i, j, k], fields.saturation_water_ice[i, j, k])
    out.decomposition_response[i, j, k] = resp
    litter_carbon = fields.litter_carbon[i, j, k]
    decomposed_litter, decomposed_fast, decomposed_slow, to_fast, to_slow, to_atmosphere =
        soil_carbon_decomposition(bgc, litter_carbon, fields.fast_carbon[i, j, k], fields.slow_carbon[i, j, k], resp)
    out.heterotrophic_respiration[i, j, k] = to_atmosphere + decomposed_fast + decomposed_slow

    # Soil-organic-matter mineralization: nitrogen released by the respired fast/slow carbon at the
    # soil C:N. The litter's own nitrogen is tracked separately (its C:N differs from the soil's).
    som_mineralization = (decomposed_fast + decomposed_slow) / bgc.soil_cn_ratio
    Δz = Δzᵃᵃᶜ(i, j, k, field_grid)
    n_release, litter_to_mineral = litter_nitrogen_exchange(
        bgc, litter_carbon, fields.litter_nitrogen[i, j, k], decomposed_litter,
        to_fast, to_slow, fields.soil_ammonium[i, j, k], Δz,
    )
    out.litter_nitrogen_release[i, j, k] = n_release
    out.net_mineralization[i, j, k] = som_mineralization + litter_to_mineral

    # Ammonia volatilization: an NH₃ sink acting only on the top soil layer (exposed to the atmosphere).
    # The LPJmL primitive takes the top-layer ammonium as a per-area amount (gN/m²) and the layer depth
    # in mm, and returns a per-day flux (gN/m²/day); convert to a per-volume per-second density rate.
    ammonium_perarea = max(zero(NF), fields.soil_ammonium[i, j, k]) * Δz * NF(1000)
    vol_flux = ammonia_volatilization(
        bgc.volatilization, fields.air_temperature[i, j], fields.windspeed[i, j], bgc.soil_ph, ammonium_perarea, Δz * NF(1000),
    )
    vol_rate = vol_flux / (NF(86400) * NF(1000) * Δz)
    out.ammonia_volatilization[i, j, k] = ifelse(k == field_grid.Nz, vol_rate, zero(NF))
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
    # Fertilizer application, distributed over the root zone like the crop fluxes.
    fertilizer_ammonium = fields.fertilizer_ammonium_flux[i, j] * per_volume
    fertilizer_nitrate = fields.fertilizer_nitrate_flux[i, j] * per_volume

    out.litter_carbon[i, j, k] = d_litter + litterfall_carbon
    out.fast_carbon[i, j, k] = d_fast
    out.slow_carbon[i, j, k] = d_slow
    # Litter nitrogen: gains the crop litterfall nitrogen (organic N, not yet mineral); loses nitrogen
    # as its carbon decomposes (routed to the soil organic matter and the mineral pool via
    # `net_mineralization`, which already carries the litter mineralization/immobilization).
    out.litter_nitrogen[i, j, k] = litterfall_nitrogen - fields.litter_nitrogen_release[i, j, k]
    # Ammonium: mineralization − immobilization (both in net_mineralization) − nitrification (both in
    # d_ammonium); crop uptake draws it down; fertilizer adds to it; NH₃ volatilizes from the top layer.
    out.soil_ammonium[i, j, k] = d_ammonium - uptake_ammonium + fertilizer_ammonium - fields.ammonia_volatilization[i, j, k]
    out.soil_nitrate[i, j, k] = d_nitrate - uptake_nitrate + fertilizer_nitrate
end
