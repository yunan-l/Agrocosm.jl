using Agrocosm
using Test

# LPJmL's crop allocation freezes `leaf` when senescence starts and only ever
# trims it through a mass-balance clamp, so a crop can finish the season holding
# a canopy's worth of carbon at zero leaf area while
# `compute_storage_carbon`'s `biomass - leaf - root` cap denies exactly that
# carbon to the grain. `senescent_leaf_release` trims `leaf` towards the
# standing canopy instead. `docs/06_allocation_bookkeeping.md` has the trace
# this was found in.

"""Build the frozen-leaf state: senescent, canopy well below the leaf pool."""
function _senescent_case(::Type{T}) where {T}
    crop = init_crop(T, 1, identity)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= Int32(1)
    crop.state.phenology.growing_days .= Int32(100)
    crop.state.phenology.senescence .= true
    crop.state.phenology.grain_set_fraction .= one(T)
    # Standing canopy is 1.0 of a potential 3.0: `lai_npp_deficit` is the
    # carbon the crop could never afford, exactly as the pre-senescence branch
    # books it.
    crop.state.canopy.lai .= T(3)
    crop.state.canopy.lai_npp_deficit .= T(2)
    crop.state.carbon.biomass .= T(200)
    # 120 gC of leaf supports LAI 120 * sla, far above the standing 1.0.
    crop.state.carbon.leaf .= T(120)
    crop.state.carbon.root .= T(20)
    crop.state.carbon.storage .= T(10)
    crop.state.carbon.pool .= T(50)
    crop.state.nitrogen.sufficiency .= one(T)
    # `compute_seasonal_nitrogen_sufficiency` reads the accumulated sum over
    # growing days as a percentage, and that percentage drives the root
    # fraction. Seed it at 100% so this case exercises the leaf release rather
    # than a 40% root allocation and the mass-balance clamp behind it.
    crop.state.nitrogen.stress_sum .= T(99)
    crop.state.water.sufficiency .= T(100)
    crop.auxiliary.phenology.fphu .= one(T)
    crop.auxiliary.stress.water_deficit .= T(100)
    # No assimilation this step, so the trace is the allocation alone.
    crop.fluxes.carbon.gross_assimilation .= zero(T)
    crop.fluxes.carbon.leaf_respiration .= zero(T)
    crop.fluxes.carbon.respiration .= zero(T)
    return crop, state
end

_release_params(::Type{T}, value) where {T} = Agrocosm.LPJmLParams{T}(;
    (f => (f === :senescent_leaf_release ? T(value) :
           getfield(Agrocosm.LPJmLParams{T}(), f))
     for f in fieldnames(Agrocosm.LPJmLParams))...)

@testset "Senescent leaf carbon release" begin
    for T in (Float32, Float64)
        sla = T(cft1.sla)

        # Zero release is LPJmL, bitwise: `leaf -= 0 * surplus` cannot perturb
        # the last bit, so this is the exact retreat the design promises.
        frozen_crop, frozen_state = _senescent_case(T)
        carbon_allocation!(cft1, frozen_state;
                           lpjmlparams = _release_params(T, 0))
        @test only(frozen_crop.state.carbon.leaf) === T(120)

        released_crop, released_state = _senescent_case(T)
        carbon_allocation!(cft1, released_state;
                           lpjmlparams = _release_params(T, 1))
        # Full release leaves exactly the standing canopy's carbon.
        standing = (T(3) - T(2)) / sla
        @test only(released_crop.state.carbon.leaf) ≈ standing

        # Half release moves exactly half the surplus in one step. The rule runs
        # every senescence day, so across a season an intermediate value decays
        # the surplus geometrically rather than releasing that share once -- the
        # per-step linearity asserted here is what makes that behaviour
        # predictable, not a claim that the parameter interpolates the outcome.
        half_crop, half_state = _senescent_case(T)
        carbon_allocation!(cft1, half_state;
                           lpjmlparams = _release_params(T, 0.5))
        @test only(half_crop.state.carbon.leaf) ≈
            T(120) - T(0.5) * (T(120) - standing)

        for crop in (frozen_crop, released_crop, half_crop)
            carbon = crop.state.carbon
            # Biomass is untouched and the four pools still close onto it: the
            # release moves carbon between pools, it does not create any.
            @test only(carbon.biomass) ≈ T(200)
            @test only(carbon.leaf) + only(carbon.root) +
                  only(carbon.storage) + only(carbon.pool) ≈ only(carbon.biomass)
            @test only(carbon.leaf) >= zero(T)
            @test only(carbon.pool) >= zero(T)
        end

        # The point of the fix: with the leaf pool frozen, storage is pinned to
        # `biomass - leaf - root` instead of to the harvest index. Releasing the
        # phantom canopy lets the harvest index take over, so storage rises
        # while biomass does not.
        frozen_store = only(frozen_crop.state.carbon.storage)
        released_store = only(released_crop.state.carbon.storage)
        @test frozen_store ≈ T(200) - T(120) - only(frozen_crop.state.carbon.root)
        @test released_store > frozen_store
        # And it stops at the harvest index rather than taking everything the
        # cap now allows.
        released_root = only(released_crop.state.carbon.root)
        @test released_store < T(200) - standing - released_root
    end
end

@testset "Senescent leaf release leaves a standing canopy alone" begin
    # A crop whose leaf pool already matches its canopy has nothing to release,
    # so the parameter must be inert there. Otherwise the fix would be silently
    # reallocating carbon in every ordinary senescence step.
    for T in (Float32, Float64)
        sla = T(cft1.sla)
        for release in (0, 1)
            crop, state = _senescent_case(T)
            crop.state.canopy.lai_npp_deficit .= zero(T)
            crop.state.carbon.leaf .= T(3) / sla
            crop.state.carbon.pool .= T(200) - T(3) / sla - T(20) - T(10)
            carbon_allocation!(cft1, state; lpjmlparams = _release_params(T, release))
            @test only(crop.state.carbon.leaf) ≈ T(3) / sla
        end
    end
end
