using Agrocosm
using Test

@testset "Crop soil ammonia volatilization" begin
    v = CropVolatilization(Float64)   # length_scale = 1.0

    # Independent recomputation of the LPJmL closed form (uncapped), then min with ammonium.
    function expected(T, wind, ph, nh4, depth; L = 1.0)
        kelvin = T + 273.15
        dissociation = 10.0^(0.05 - 2788 / kelvin)
        aqueous_fraction = 1 / (1 + 10.0^(-ph) / dissociation)
        aqueous_nh3 = aqueous_fraction * nh4 / depth * 1000
        henry = 0.2138 / kelvin * 10.0^(6.123 - 1825 / kelvin)
        mt = 0.000612 * wind^0.8 * kelvin^0.382 * L^(-0.2)
        return clamp(86400 * mt * henry * aqueous_nh3, 0.0, nh4)
    end

    @testset "closed form and cap" begin
        f = Agrocosm.ammonia_volatilization(v, 20.0, 2.0, 7.0, 0.01, 0.2)
        @test f ≈ expected(20.0, 2.0, 7.0, 0.01, 0.2) rtol = 1e-9
        @test 0.0 ≤ f ≤ 0.01                         # capped at the ammonium stock
    end

    @testset "no ammonium, no flux" begin
        @test Agrocosm.ammonia_volatilization(v, 20.0, 2.0, 7.0, 0.0, 0.2) == 0.0
    end

    @testset "higher pH mobilizes more NH3 (before the cap)" begin
        # Small ammonium so the flux is not fully capped; higher pH → larger aqueous NH3 fraction.
        nh4 = 1.0e-6
        f_low = Agrocosm.ammonia_volatilization(v, 20.0, 2.0, 6.0, nh4, 0.2)
        f_high = Agrocosm.ammonia_volatilization(v, 20.0, 2.0, 8.0, nh4, 0.2)
        @test f_high ≥ f_low
    end
end
