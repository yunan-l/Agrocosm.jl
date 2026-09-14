using Test
using Agrocosm
using Agrocosm: uptake_weight, compute_layer_transpiration, cft3

@testset "the uptake weight is LPJmL until it is asked otherwise" begin
    # Identity, not approximate equality: the ablation contract needs the
    # shipped exponent to return the same bits.
    for water in (0.0, 0.079, 0.5, 0.714, 1.0)
        @test uptake_weight(water, 1.0) === water
    end
    @test cft3.uptake_availability_exponent == 1.0
    # And the layer routine with the default argument is the old expression.
    @test compute_layer_transpiration(3.0, 0.32, 0.079, 48.2)[1] ===
          3.0 * 0.32 * 0.079
end

@testset "a lower exponent stops a drying layer from being abandoned" begin
    # The measured failure: 0-20 cm at 0.714 of available water supplied 207.7 mm
    # of the season's 236, while 20-50 cm at 0.079 supplied 11.3. The weighting
    # is what concentrates it, so flattening it must move uptake downward.
    roots = (0.384, 0.319, 0.209)
    water = (0.714, 0.079, 0.831)
    share(exponent) = begin
        weights = ntuple(i -> roots[i] * uptake_weight(water[i], exponent), 3)
        total = sum(weights)
        ntuple(i -> weights[i] / total, 3)
    end
    lpjml, flat = share(1.0), share(0.0)
    # The middle layer holds 32% of the roots and supplies 5% of the water.
    @test lpjml[2] < roots[2] / 5
    @test flat[2] > 3 * lpjml[2]    # and comes back when the weighting flattens
    @test sum(lpjml) ≈ 1.0
    @test sum(flat) ≈ 1.0
    # Root density alone, which is what exponent 0 must mean.
    @test flat[1] ≈ roots[1] / sum(roots)
end

@testset "the cap is the layer's real water whatever the weighting" begin
    # A layer holding 1 mm cannot give 3 mm however the weights are set: the
    # exponent decides WHERE water is drawn from, never how much a layer holds.
    for exponent in (1.0, 0.5, 0.0)
        drawn, capped = compute_layer_transpiration(1000.0, 0.32, 0.02, 48.2, exponent)
        @test capped
        @test drawn ≈ 0.02 * 48.2
    end
end
