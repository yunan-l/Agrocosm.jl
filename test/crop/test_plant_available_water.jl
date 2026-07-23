using Agrocosm
using Test

@testset "Crop plant-available-water stress" begin
    β(θ) = Agrocosm.soil_moisture_limiting_factor(θ, 0.1, 0.4)   # wilting 0.1, field capacity 0.4

    @testset "soil-moisture limiting factor" begin
        @test β(0.1) ≈ 0.0            # at wilting point
        @test β(0.4) ≈ 1.0            # at field capacity
        @test β(0.25) ≈ 0.5           # midpoint
        @test β(0.05) == 0.0          # below wilting → clamped
        @test β(0.6) == 1.0           # above field capacity → clamped
        vals = [β(θ) for θ in 0.1:0.03:0.4]
        @test issorted(vals)          # monotone in water content
    end

    @testset "plant-available water in a layer" begin
        # (θ − θ_wilting)·Δz, floored at zero.
        @test Agrocosm.plant_available_water(0.3, 0.1, 0.5) ≈ 0.1
        @test Agrocosm.plant_available_water(0.05, 0.1, 0.5) == 0.0   # below wilting
    end
end
