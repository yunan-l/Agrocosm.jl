using Agrocosm
using Test

@testset "Crop nitrogen demand" begin
    d = CropNitrogenDemand(Float64)
    nd(v, lc, rc, pc, sc, T) = Agrocosm.crop_nitrogen_demand(d, v, lc, rc, pc, sc, T)

    @testset "closed form at 25 °C" begin
        dl, dt = nd(10.0, 100.0, 50.0, 10.0, 5.0, 25.0)
        rubisco = 25.0 * 1.0e-3 * 10.0 / (86400 * 12 * 1.0e-6)   # exp(0) at 25 °C
        demand_leaf = rubisco + (1 / 58.8) * 100.0
        @test dl ≈ demand_leaf rtol = 1e-9
        nc = clamp(demand_leaf / 100.0, 1 / 58.8, 1 / 14.3)
        demand_total = demand_leaf + nc * (50.0 / 1.16 + 10.0 / 3.0 + 5.0 / 0.99)
        @test dt ≈ demand_total rtol = 1e-9
        @test dt > dl                       # total exceeds leaf-only demand
    end

    @testset "demand inverts the Vcmax capacity" begin
        # The leaf Rubisco demand for a given Vcmax equals the Vcmax that that much nitrogen supports:
        # feeding demand_leaf back through the nitrogen limitation recovers the Vcmax (leaf N above
        # structural minimum). Uses matching parameters.
        lim = CropNitrogenVcmaxLimit(Float64)
        dl, _ = nd(10.0, 100.0, 0.0, 0.0, 0.0, 25.0)
        v, _ = Agrocosm.nitrogen_limited_vcmax(lim, 1000.0, dl, 100.0, 25.0)  # ample potential
        @test v ≈ 10.0 rtol = 1e-6
    end

    @testset "temperature and carbon dependence" begin
        dl_cold, _ = nd(10.0, 100.0, 50.0, 10.0, 5.0, 15.0)
        dl_warm, _ = nd(10.0, 100.0, 50.0, 10.0, 5.0, 35.0)
        @test dl_warm < dl_cold             # warmer → lower Rubisco N requirement
        # zero leaf carbon: nc_ratio clamps to the minimum, leaf demand is Rubisco-only
        dl0, dt0 = nd(10.0, 0.0, 50.0, 0.0, 0.0, 25.0)
        rubisco = 25.0 * 1.0e-3 * 10.0 / (86400 * 12 * 1.0e-6)
        @test dl0 ≈ rubisco rtol = 1e-9
        @test dt0 ≈ rubisco + (1 / 58.8) * (50.0 / 1.16) rtol = 1e-9
    end
end
