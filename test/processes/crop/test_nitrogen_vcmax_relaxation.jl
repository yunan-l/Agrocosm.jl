using Agrocosm
using Test

@testset "Lifting the nitrogen clamp is isolated and zero is the clamp bitwise" begin
    function run_clamp(relaxation; leaf_nitrogen = 0.05f0, leaf_carbon = 20.0f0,
                       potential = 40.0f0)
        crop = init_crop(1, identity)
        crop.state.phenology.is_growing .= 1
        crop.auxiliary.photosynthesis.potential_vcmax .= potential
        crop.auxiliary.photosynthesis.lambda .= 0.7f0
        crop.auxiliary.stress.nitrogen_demand_leaf .= leaf_nitrogen
        crop.state.carbon.leaf .= leaf_carbon
        state = test_model_state(crop, init_soil(1, soilparams.soildepth, identity))
        params = Agrocosm.LPJmLParams{Float32}(;
            (f => (f === :nitrogen_vcmax_relaxation ? Float32(relaxation) :
                   getfield(Agrocosm.LPJmLParams{Float32}(), f))
             for f in fieldnames(Agrocosm.LPJmLParams))...)
        limit_vcmax_by_nitrogen!(state, cft1, Float32[20]; lpjmlparams = params)
        return (vcmax = crop.auxiliary.photosynthesis.vcmax[1],
                ratio = crop.auxiliary.photosynthesis.nitrogen_limitation[1])
    end

    # CONTRACT. Zero must reproduce the shipped clamp exactly, not approximately.
    shipped = let crop = init_crop(1, identity)
        crop.state.phenology.is_growing .= 1
        crop.auxiliary.photosynthesis.potential_vcmax .= 40.0f0
        crop.auxiliary.photosynthesis.lambda .= 0.7f0
        crop.auxiliary.stress.nitrogen_demand_leaf .= 0.05f0
        crop.state.carbon.leaf .= 20.0f0
        state = test_model_state(crop, init_soil(1, soilparams.soildepth, identity))
        limit_vcmax_by_nitrogen!(state, cft1, Float32[20])
        crop.auxiliary.photosynthesis.vcmax[1]
    end
    @test run_clamp(0.0).vcmax === shipped
    @test shipped < 40.0f0        # the fixture must actually be clamped, or the test is vacuous

    # Full relaxation removes the clamp from photosynthesis entirely.
    @test run_clamp(1.0).vcmax ≈ 40.0f0
    @test run_clamp(1.0).ratio ≈ 1.0f0

    # Monotone in between, and the half-way point is the midpoint of the blend.
    previous = shipped
    for relaxation in (0.0f0, 0.25f0, 0.5f0, 0.75f0, 1.0f0)
        current = run_clamp(relaxation).vcmax
        @test current >= previous
        previous = current
    end
    @test run_clamp(0.5).vcmax ≈ shipped + 0.5f0 * (40.0f0 - shipped)

    # Out-of-range values are clamped rather than extrapolating past `potential`.
    @test run_clamp(3.0).vcmax ≈ run_clamp(1.0).vcmax
    @test run_clamp(-1.0).vcmax === shipped

    # A crop with ample leaf nitrogen is not clamped, so relaxation changes
    # nothing there - the parameter must not be a blanket vcmax multiplier.
    ample = run_clamp(0.0; leaf_nitrogen = 50.0f0)
    @test ample.vcmax ≈ 40.0f0
    @test run_clamp(1.0; leaf_nitrogen = 50.0f0).vcmax ≈ ample.vcmax
end
