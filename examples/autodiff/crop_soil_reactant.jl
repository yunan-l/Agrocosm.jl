# # Running the crop soil biogeochemistry on Reactant
#
# Terrarium (and therefore Agrocosm) can trace and compile a model run through
# [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl) → MLIR/StableHLO → XLA, which is how the
# model runs fast on GPUs and, later, how the whole rollout is differentiated in one compiled program
# (see `differentiating_crop_soil_reactant.jl`). The only user-facing change is the architecture: we
# build the grid on `ReactantState()` and the state lives on the device; `run!` traces the stepping
# loop and compiles it once.
#
# Run from this environment:  julia --project=examples/autodiff examples/autodiff/crop_soil_reactant.jl

using Agrocosm, Terrarium
using Reactant, CUDA   # CUDA provides Reactant's kernel integration, even on CPU

NF = Float32

# Identical model to the CPU example — only the architecture differs.
grid = ColumnGrid(ReactantState(), NF, UniformSpacing(Δz = NF(0.1), N = 10))
soil = SoilEnergyWaterCarbon(NF; biogeochem = CropSoilBiogeochemistry(NF))
model = SoilModel(grid; soil)
bcs = PrescribedSurfaceTemperature(:T_ub, NF(15))
initializers = (temperature = (x, z) -> NF(15) - NF(0.02) * z,)
integrator = initialize(model; boundary_conditions = bcs, initializers)

# `run!` on a `ReactantState` integrator compiles the stepping loop the first time it is called (a
# few seconds to minutes) and then executes the compiled program. A week of hourly-ish steps lets the
# soil carbon–nitrogen pools evolve visibly from their initial state.
Terrarium.run!(integrator; steps = 1000, Δt = NF(600))

# Materialize a few soil-column diagnostics back onto the host (`on_architecture(CPU(), …)` is the
# canonical transfer). The crop C–N pools have evolved from their initial state.
using Oceananigans: on_architecture
using Oceananigans.Architectures: CPU
host(field) = Array(interior(on_architecture(CPU(), field)))[1, 1, :]

println("soil nitrate  (kgN/m³) by depth: ", round.(host(integrator.state.soil_nitrate), sigdigits = 3))
println("soil ammonium (kgN/m³) by depth: ", round.(host(integrator.state.soil_ammonium), sigdigits = 3))
println("litter carbon (kgC/m³) by depth: ", round.(host(integrator.state.litter_carbon), sigdigits = 3))
