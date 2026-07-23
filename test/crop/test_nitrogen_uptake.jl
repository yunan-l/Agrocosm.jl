using Agrocosm
using Test

@testset "Crop root nitrogen uptake kinetics" begin
    @testset "temperature response" begin
        f(T) = Agrocosm.nitrogen_uptake_temperature_response(T, -25.0, 15.0, 15.0)
        @test f(15.0) ≈ 1.0            # normalized to 1 at the reference temperature
        @test f(-25.0) ≈ 0.0           # zero at T_0
        @test f(-30.0) == 0.0          # below T_0 → clamped to 0
        @test f(15.0) > f(0.0) > f(-20.0)   # rises from T_0 toward the optimum
    end

    @testset "Michaelis-Menten uptake potential" begin
        k = CropNitrogenUptakeKinetics(Float64)   # vmax=1.5, kmin=0.05, Km=0.70
        up(N; scale = 1.0, rf = 1.0) = Agrocosm.root_nitrogen_uptake_potential(k, N, scale, rf)

        @test up(0.0) == 0.0                       # empty pool → no uptake
        # saturation = N/(N + Km·scale); at N = Km·scale it is 0.5.
        # potential = vmax·(kmin + saturation)·rf, capped at the available N. Use a small root
        # factor so the cap does not bind and the kinetic form is exercised directly.
        N = 0.70                                   # = Km·scale with scale=1 → saturation 0.5
        pot = 1.5 * (0.05 + 0.5) * 0.5             # rf = 0.5 → potential = 0.4125 < N
        @test Agrocosm.root_nitrogen_uptake_potential(k, N, 1.0, 0.5) ≈ pot rtol = 1e-9
        # monotone increasing in available nitrogen via saturation (small rf, cap not binding)
        vals = [up(N; rf = 0.1) for N in 0.1:0.2:2.0]
        @test issorted(vals)
        # capped at available nitrogen when the root factor is large
        @test up(0.5; rf = 1.0e6) ≈ 0.5
        # scales linearly with the root factor when the cap does not bind
        @test up(1.0; rf = 0.6) ≈ 2 * up(1.0; rf = 0.3) rtol = 1e-12
    end
end
