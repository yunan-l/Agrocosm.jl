# Minimal pure-Terrarium reproduction of the Reactant incompatibility that blocks compiling the full
# crop `LandModel` (see 2026-07-24_NOTES_future_work.md → "Reactant compatibility of the crop
# LandModel"). It uses NO Agrocosm code — only Terrarium's own `StaticExponentialRootDistribution` —
# so it can be dropped into Terrarium's `test/reactant/` when fixing this upstream.
#
# Symptom: materializing the static root-fraction field on a `ReactantState` grid throws
#   MethodError: no method matching exp(::Reactant.TracedRArray{Float32, 3})
# The normalization `sum(R, dims = 3)` over the `∂R∂z` `FunctionField` (whose closure calls the scalar
# `exp` in `root_density`) evaluates that closure on the whole traced z-array instead of per element, so
# `exp` is applied to a `TracedRArray`. Agrocosm's `CropRootDistribution.crop_root_fraction` mirrors the
# same Terrarium pattern (`FunctionField` × `Δz`, `/ sum(·, dims = 3)`), so it fails the same way, which
# is why the crop `LandModel` cannot yet be built/compiled on `ReactantState`.
#
# Terrarium's Reactant suite currently exercises only `SoilModel` configs (no vegetation), so this path
# is untested upstream.
#
# Run (Julia 1.12, Reactant-enabled env):  julia +1.12 --project=<env> reactant_rootfraction_repro.jl

using Terrarium
using Reactant, CUDA

NF = Float32
grid = ColumnGrid(ReactantState(), NF, ExponentialSpacing(Δz_min = 0.05f0, Δz_max = 1.0f0, N = 20))
rootdist = Terrarium.StaticExponentialRootDistribution(NF)

println("Materializing StaticExponentialRootDistribution root_fraction on ReactantState...")
try
    R_norm = Terrarium.root_fraction(rootdist, grid, nothing, (;))
    field = Field(R_norm)
    compute!(field)
    println("RESULT: OK — root_fraction materialized: ", Array(interior(field))[1, 1, 1:3])
catch e
    io = IOBuffer()
    showerror(io, e)
    println("RESULT: FAILED — ", first(String(take!(io)), 400))
end
