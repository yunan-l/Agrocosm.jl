using Test
using Agrocosm

# The two diagnostics docs/16 turns on. Neither may touch the physics, and the
# harvest-index one has to mean exactly what its name says, or the measurement it
# exists to make is worthless.

@testset "harvest_index_binds agrees with compute_storage_carbon" begin
    # `harvest_index_binds` duplicates `compute_storage_carbon`'s expressions
    # rather than refactoring them out, so that function stays bitwise untouched
    # on the differentiated path. Duplication drifts; this pins it.
    #
    # The property under test is the flag's MEANING: it is true exactly when
    # moving the harvest index moves storage carbon. That is what decides whether
    # a mechanism multiplying the index can matter on a given day, which is the
    # whole question docs/16 asks.
    T = Float64
    storage(biomass, leaf, root, froot, hi, hiopt, deposited) =
        max(deposited, Agrocosm.compute_storage_carbon(
            T(biomass), T(leaf), T(root), T(froot), T(hi), T(hiopt),
        ))

    rng_state = 12345
    next() = (rng_state = (1103515245 * rng_state + 12345) % 2147483648;
              rng_state / 2147483648)

    # The probe reduces the harvest index by a fixed relative step. That step can
    # itself carry `candidate` across `mass_cap` on a draw where it started just
    # above it: the flag correctly says the index did not bind, yet storage moves
    # because the reduced index made it bind. Those draws are counted and
    # required to BE boundary crossings - they are a property of the finite
    # probe, not of the predicate, and silently widening the tolerance would hide
    # a real drift just as effectively.
    step = 0.999
    agreements = 0
    binding = 0
    boundary = 0
    for _ in 1:20000
        biomass = 10 * next()
        leaf = 5 * next()
        root = 5 * next()
        froot = 0.9 * next()
        hi = 0.05 + 0.9 * next()
        # Both branches of `optimal_index > one(T)`: cft6 (hiopt 3.5) and cft7
        # (hiopt 2) take the second one, every other CFT the first.
        hiopt = next() < 0.25 ? 1 + 3 * next() : 0.2 + 0.7 * next()
        deposited = 3 * next()

        flag = Agrocosm.harvest_index_binds(
            T(biomass), T(leaf), T(root), T(froot), T(hi), T(hiopt), T(deposited),
        )
        nudged_flag = Agrocosm.harvest_index_binds(
            T(biomass), T(leaf), T(root), T(froot), T(hi * step), T(hiopt), T(deposited),
        )
        base = storage(biomass, leaf, root, froot, hi, hiopt, deposited)
        nudged = storage(biomass, leaf, root, froot, hi * step, hiopt, deposited)
        moved = base != nudged

        if flag == moved
            agreements += 1
        elseif flag != nudged_flag
            boundary += 1          # the probe itself changed which term bound
        end
        flag && (binding += 1)
    end
    @test agreements + boundary == 20000
    @test boundary < 20            # a rare edge, not the common case
    # A test where the flag is never true would pass vacuously.
    @test binding > 2000
end

@testset "the flowering window gates window NPP, and does not weight it" begin
    # The damage mechanisms weight the window with a raised cosine. This
    # diagnostic must NOT: grain number responds to assimilate supply over the
    # critical period, so a day just inside the window contributes its whole NPP.
    T = Float32
    start, stop = T(0.45), T(0.70)
    inside(fphu) = Agrocosm.flowering_weight(T(fphu), start, stop) > zero(T)

    @test !inside(0.20)          # before the window
    @test !inside(0.90)          # after it
    @test !inside(0.45)          # open at the lower end
    @test !inside(0.70)          # open at the upper end
    @test inside(0.46)
    @test inside(0.575)
    @test inside(0.69)
    # A degenerate window admits nothing rather than everything.
    @test Agrocosm.flowering_weight(T(0.5), T(0.6), T(0.6)) == zero(T)
    @test Agrocosm.flowering_weight(T(0.5), T(0.7), T(0.6)) == zero(T)
    # The gate is a threshold on the weight, so the weight's shape inside the
    # window must not reach it: the smallest interior weight is still positive.
    @test Agrocosm.flowering_weight(T(0.4501), start, stop) > zero(T)
end

@testset "S1/S2: the diagnostics are inert on the physics" begin
    # Neither diagnostic feeds any model state. The contract is stated at
    # `accumulate_season_process_diagnostics!` and rests on both arrays being
    # written only by kernels nothing else reads - if that ever stops being true,
    # every ablation rung's bitwise-equivalence claim goes with it.
    crop = init_crop(4, identity)
    stress = crop.auxiliary.stress
    @test hasproperty(stress, :harvest_index_binding)
    @test all(iszero, stress.harvest_index_binding)

    output = init_output(4, identity)
    for field in (:window_npp, :hi_binding_days)
        @test hasproperty(output.crop, field)
        @test hasproperty(output.annual, field)   # harvesting.jl indexes both by one symbol
    end
    for field in (:active_window_npp, :active_hi_binding_days)
        @test hasproperty(output.annual, field)
        @test all(iszero, getproperty(output.annual, field))
    end
end
