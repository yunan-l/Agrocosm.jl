# Crop soil carbon pool dynamics (LPJmL `soil_carbon`). Each pool (litter, fast, slow) decomposes by
# a first-order rate scaled by the environmental decomposition response; decomposed litter carbon is
# partitioned into the fast and slow soil pools and the atmosphere (heterotrophic respiration), and
# decomposed fast/slow carbon is respired. Tested scalar physics; the full multi-layer pool state and
# its coupling to the soil nitrogen pools are assembled in the crop soil biogeochemistry (plan Phase 5).

"""
    $(TYPEDEF)

Crop soil carbon decomposition parameters (daily first-order rates at 10 °C and litter routing
fractions).

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropSoilCarbon{NF}
    "Fast-pool decomposition rate at 10 °C"
    @param k_fast::NF = 0.04 / 365 (units = u"d^-1", bounds = Positive)
    "Slow-pool decomposition rate at 10 °C"
    @param k_slow::NF = 0.001 / 365 (units = u"d^-1", bounds = Positive)
    "Fraction of retained decomposed litter routed to the fast pool"
    @param fast_fraction::NF = 0.98 (bounds = UnitInterval,)
    "Fraction of decomposed litter emitted directly to the atmosphere"
    @param atmospheric_fraction::NF = 0.5 (bounds = UnitInterval,)
end

CropSoilCarbon(::Type{NF}; kwargs...) where {NF} = CropSoilCarbon{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

First-order decomposed carbon `(1 − exp(−rate·response))·pool`, the fraction of a carbon pool
decomposed in one step at the given decomposition `response`. Uses `-expm1` for accuracy near zero.
"""
@inline decomposed_carbon(rate::NF, response::NF, pool::NF) where {NF} = -expm1(-rate * response) * pool

"""
    $(TYPEDSIGNATURES)

Route decomposed litter carbon into `(to_fast, to_slow, to_atmosphere)`. The retained fraction
`(1 − atmospheric_fraction)` is split between fast and slow by `fast_fraction`; the rest is respired.
The three components sum to `decomposed_litter`.
"""
@inline function route_litter_carbon(c::CropSoilCarbon{NF}, decomposed_litter::NF) where {NF}
    retained = (one(NF) - c.atmospheric_fraction) * decomposed_litter
    to_fast = c.fast_fraction * retained
    to_slow = (one(NF) - c.fast_fraction) * retained
    to_atmosphere = c.atmospheric_fraction * decomposed_litter
    return to_fast, to_slow, to_atmosphere
end

"""
    $(TYPEDSIGNATURES)

Heterotrophic respiration from decomposing the three soil carbon pools: the atmospheric fraction of
the decomposed litter plus all decomposed fast and slow carbon.
"""
@inline function heterotrophic_respiration(
        c::CropSoilCarbon{NF}, decomposed_litter::NF, decomposed_fast::NF, decomposed_slow::NF,
    ) where {NF}
    return c.atmospheric_fraction * decomposed_litter + decomposed_fast + decomposed_slow
end
