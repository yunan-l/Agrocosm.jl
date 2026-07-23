"""
    $(TYPEDEF)

LPJmL-style crop root distribution. The cumulative fraction of roots above depth
`d` (measured downward from the surface) is

```math
Y(d) = 1 - β^{d}
```

where `β = beta_root` (dimensionless, `0 < β < 1`) and `d` is depth in centimetres,
following LPJmL (`rootdist` in `initgrid.c`/`newpft.c`). The corresponding
continuous root density with respect to the vertical coordinate `z` (metres,
negative downward, so `d = -100 z`) is

```math
\\frac{\\partial R}{\\partial z} \\propto β^{d} = \\exp(-100\\, z \\ln β),
```

which Terrarium integrates over the soil column (midpoint rule) and normalizes to
sum to unity. Because normalization removes any constant prefactor, the discretized
per-layer fractions reproduce the LPJmL layer fractions
`(β^{d_{l-1}} - β^{d_l}) / (1 - β^{d_{bottom}})` in the fine-grid limit.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropRootDistribution{NF} <: Terrarium.AbstractRootDistribution{NF}
    "LPJmL exponential root-profile parameter β (per-centimetre base, `0 < β < 1`)"
    @param beta_root::NF = 0.94 (bounds = UnitInterval,)
end

CropRootDistribution(::Type{NF}; kwargs...) where {NF} = CropRootDistribution{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Continuous root density `∂R/∂z` at height `z` (metres, negative downward). The
result carries an arbitrary positive constant factor that is removed by the
column normalization in [`crop_root_fraction`](@ref).
"""
@inline function Terrarium.root_density(rd::CropRootDistribution{NF}, z) where {NF}
    lnβ = log(rd.beta_root)
    depth_cm = -NF(100) * z
    return -lnβ * exp(depth_cm * lnβ)
end

Terrarium.variables(rd::CropRootDistribution) = (
    Terrarium.auxiliary(:root_fraction, XYZ(), crop_root_fraction, rd),
)

"""
    $(TYPEDSIGNATURES)

Return a normalized root-fraction field over a column grid: the continuous
[`Terrarium.root_density`](@ref) sampled at cell centres, weighted by the layer
thickness, and normalized to sum to unity over the column.
"""
function crop_root_fraction(rd::CropRootDistribution{NF}, grid::Terrarium.AbstractColumnGrid, clock, fields) where {NF}
    fgrid = get_field_grid(grid)
    ∂R∂z = FunctionField{Center, Center, Center}(fgrid, parameters = rd) do x, z, params
        Terrarium.root_density(params, z)
    end
    Δz = zspacings(fgrid, Center(), Center(), Center())
    R = ∂R∂z * Δz
    R_norm = R / sum(R, dims = 3)
    return R_norm
end
