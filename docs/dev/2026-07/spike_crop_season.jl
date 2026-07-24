# Phase 5 seasonal spike: run the crop VegetationModel through a full phenological cycle and check
# that the leaf area index rises to a peak near the senescence onset (fphu ≈ fphusen) and then
# declines to near zero at maturity (fphu → 1) — the LPJmL LAI trajectory, produced dynamically by
# the prognostic heat-unit accumulation.
#
# Run: julia --project=. docs/dev/2026-07/spike_crop_season.jl

using Agrocosm
using Terrarium

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
soil = SoilEnergyWaterCarbon(eltype(grid))

# Small heat-unit requirement so the crop completes its cycle in a short spike.
vegetation = CropVegetation(
    eltype(grid);
    phenology_dynamics = CropPhenologyDynamics(eltype(grid); heat_unit_requirement = 60.0),
)
land = LandModel(grid; soil, vegetation)

function run_season(land)
    integrator = initialize(land; initializers = (temperature = 20.0,))
    set!(integrator.state.air_temperature, 25.0)
    set!(integrator.state.surface_shortwave_down, 400.0)

    Δt = 600.0
    peak_lai = 0.0
    peak_fphu = 0.0
    trajectory = Tuple{Float64, Float64}[]
    for _ in 1:12
        run!(integrator; steps = 30, Δt = Δt)
        fphu = interior(integrator.state.phenology_heat_unit_fraction)[1, 1, 1]
        lai = interior(integrator.state.leaf_area_index)[1, 1, 1]
        push!(trajectory, (fphu, lai))
        if lai > peak_lai
            peak_lai = lai
            peak_fphu = fphu
        end
    end
    return trajectory, peak_lai, peak_fphu
end

trajectory, peak_lai, peak_fphu = run_season(land)
final_fphu, final_lai = trajectory[end]

println("SPIKE OK")
println("  LAI trajectory (fphu, LAI):")
for (fphu, lai) in trajectory
    println("    fphu=", round(fphu, digits = 3), "  LAI=", round(lai, digits = 3))
end
println("  peak LAI = ", round(peak_lai, digits = 3), " at fphu = ", round(peak_fphu, digits = 3))
println("  final: fphu = ", round(final_fphu, digits = 3), ", LAI = ", round(final_lai, digits = 3))

@assert peak_lai > 1.0 "the canopy should develop a substantial peak LAI"
@assert 0.5 < peak_fphu < 0.85 "peak LAI should occur near the senescence onset (fphusen = 0.7)"
@assert final_fphu > 0.9 "the crop should approach maturity"
@assert final_lai < peak_lai "LAI should decline from its peak during senescence"
println("SPIKE ASSERTIONS PASSED")
