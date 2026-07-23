using Agrocosm
using Test

@testset "Crop soil nitrification" begin
    n = CropNitrification(Float64)   # k_max=0.10, k_2=0.01, a=0.45,b=1.27,c=0.0012,d=2.84

    @testset "temperature factor (Gaussian, peak 18.79 °C)" begin
        @test Agrocosm.nitrification_temperature_factor(n, 18.79) ≈ 1.0
        @test Agrocosm.nitrification_temperature_factor(n, 5.0) < 1.0
        @test Agrocosm.nitrification_temperature_factor(n, 18.79) > Agrocosm.nitrification_temperature_factor(n, 40.0)
    end

    @testset "pH factor" begin
        @test Agrocosm.nitrification_ph_factor(n, 5.0) ≈ 0.56       # atan(0)=0 at pH 5
        # monotone increasing in pH
        @test Agrocosm.nitrification_ph_factor(n, 7.0) > Agrocosm.nitrification_ph_factor(n, 4.0)
    end

    @testset "moisture factor (peaked, closed form)" begin
        # matches the LPJmL closed form on its support
        wfps = 0.6
        n_nit = 0.45 - 1.27
        m_nit = 0.45 - 0.0012
        z_nit = 2.84 * (1.27 - 0.45) / (0.45 - 0.0012)
        b1 = (wfps - 1.27) / n_nit
        b2 = (wfps - 0.0012) / m_nit
        @test Agrocosm.nitrification_moisture_factor(n, wfps) ≈ b1^z_nit * b2^2.84 rtol = 1e-9
        @test Agrocosm.nitrification_moisture_factor(n, 0.0) ≥ 0.0
        @test Agrocosm.nitrification_moisture_factor(n, 1.0) ≥ 0.0
    end

    @testset "gross nitrification capped by ammonium, N2O split" begin
        g, n2o = Agrocosm.gross_nitrification(n, 100.0, 0.6, 18.79, 7.0)
        @test 0.0 ≤ g ≤ 100.0
        @test n2o ≈ 0.01 * g
        # proportional to ammonium at the same environment
        g2, _ = Agrocosm.gross_nitrification(n, 50.0, 0.6, 18.79, 7.0)
        @test g2 ≈ g / 2 rtol = 1e-9
        # cannot exceed the ammonium stock even with maximal factors
        gcap, _ = Agrocosm.gross_nitrification(n, 1.0, 0.6, 18.79, 14.0)
        @test gcap ≤ 1.0
    end
end
