# Model registry for the CPU-vs-Reactant correctness suite (mirrors Terrarium's test/reactant/setup.jl).
#
# Add a tested configuration with a `build_model(::Val{:name}, arch, NF)` method returning a NamedTuple
# `(; model, boundary_conditions, initializers, Δt)`. `arch` is the ONLY thing that differs between the
# CPU and Reactant runs — everything else is identical.

using Agrocosm

# --- generic helpers --------------------------------------------------------------------

function build_integrator(v::Val, arch, NF)
    cfg = build_model(v, arch, NF)
    return Terrarium.initialize(
        cfg.model;
        boundary_conditions = cfg.boundary_conditions,
        initializers = cfg.initializers,
    )
end

cpu_dt(v::Val, NF) = build_model(v, CPU(), NF).Δt

# --- :crop_soil_biogeochemistry — soil column carrying the crop C–N biogeochemistry -----
# A `SoilModel` whose `biogeochem` slot is Agrocosm's `CropSoilBiogeochemistry` (prognostic
# litter/fast/slow carbon + ammonium/nitrate, with mineralization/nitrification/denitrification). The
# surface is a warm prescribed temperature so decomposition and the mineral-N transforms are active.

function build_model(::Val{:crop_soil_biogeochemistry}, arch, NF)
    grid = ColumnGrid(arch, NF, UniformSpacing(Δz = NF(0.1), N = 10))
    soil = SoilEnergyWaterCarbon(NF; biogeochem = CropSoilBiogeochemistry(NF))
    model = SoilModel(grid; soil)
    bcs = PrescribedSurfaceTemperature(:T_ub, NF(15))
    inits = (temperature = (x, z) -> NF(15) - NF(0.02) * z,)
    return (; model, boundary_conditions = bcs, initializers = inits, Δt = NF(600))
end

# --- :crop_soil_biogeochemistry_stretched — same, on an ExponentialSpacing grid ---------
# Array-valued (stretched) vertical coordinates exercise the array-z tracing path under Reactant.

function build_model(::Val{:crop_soil_biogeochemistry_stretched}, arch, NF)
    grid = ColumnGrid(arch, NF, ExponentialSpacing(Δz_min = NF(0.05), Δz_max = NF(1), N = 10))
    soil = SoilEnergyWaterCarbon(NF; biogeochem = CropSoilBiogeochemistry(NF))
    model = SoilModel(grid; soil)
    bcs = PrescribedSurfaceTemperature(:T_ub, NF(15))
    inits = (temperature = (x, z) -> NF(15) - NF(0.02) * z,)
    return (; model, boundary_conditions = bcs, initializers = inits, Δt = NF(600))
end
