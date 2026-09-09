using Agrocosm
using Test

@testset "guarded_quotient" begin
    for T in (Float32, Float64)
        # Bitwise identical to `/` wherever the denominator is at least `eps`.
        for (n, d) in ((T(3), T(2)), (T(-7), T(11)), (T(0.5), T(1e-3)),
                       (T(1), eps(T)), (T(0), T(4)))
            @test guarded_quotient(n, d) === n / d
        end
        # Finite at exactly zero, where `/` is not. The value is not meaningful
        # there -- callers keep their own `> 0` guard for that -- only finite.
        @test isfinite(guarded_quotient(one(T), zero(T)))
        @test guarded_quotient(zero(T), zero(T)) === zero(T)
        @test guarded_quotient(one(T), zero(T)) === one(T) / eps(T)
        # Monotone and sign-preserving, so a caller's clamp still behaves.
        @test guarded_quotient(-one(T), zero(T)) < zero(T)
        @test guarded_quotient(one(T), T(1e-3)) > guarded_quotient(one(T), T(1e-2))
    end
end
