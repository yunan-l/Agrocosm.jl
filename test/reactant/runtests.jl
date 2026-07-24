# CPU-vs-Reactant correctness + reverse-mode autodiff tests for Agrocosm's crop processes (mirrors
# Terrarium's test/reactant/runtests.jl).
#
# Reactant is NOT a dependency of the main test suite. Run from the repo root with this dedicated
# environment, on Julia 1.12:
#
#     julia +1.12 --project=test/reactant test/reactant/runtests.jl
#
# Reactant raising requires @inbounds elision: with --check-bounds=yes (as Pkg.test forces by default)
# every kernel retains its bounds-check throw paths, which lower to llvm.intr.trap in the GPU pipeline
# and cannot be raised to StableHLO — every compile would fail with an opaque MLIR error. Fail fast
# with a readable message instead (checked before any `using`, so we do not first recompile the stack).
if Base.JLOptions().check_bounds == 1  # 0 = auto (default), 1 = yes, 2 = no
    error(
        "The Reactant correctness tests cannot run with forced bounds checking " *
            "(--check-bounds=yes): bounds checks make every kernel un-raisable under " *
            "Reactant. Run directly with default flags: " *
            "julia +1.12 --project=test/reactant test/reactant/runtests.jl"
    )
end

using Agrocosm
using Terrarium
using Reactant
using CUDA   # required by Reactant's KernelAbstractions integration, even on CPU
using Test

# Tolerances: XLA reorders floating-point ops, so exact equality is not expected.
const DEFAULT_NF = Float32
const NSTEPS = 100
const RTOL = 1.0e-3
const ATOL = 1.0e-6

include("correctness.jl")
include("setup.jl")

@testset "Agrocosm CPU vs Reactant" begin
    test_model(:crop_soil_biogeochemistry)
    test_model(:crop_soil_biogeochemistry_stretched)   # array-valued (ExponentialSpacing) vertical coordinates
end

include("autodiff.jl")
