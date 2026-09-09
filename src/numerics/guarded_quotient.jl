"""
    guarded_quotient(numerator, denominator)

Divide with the denominator floored at `eps(T)`, so the quotient and **its
reverse-mode tangent** stay finite when the denominator reaches exactly zero.

This exists because the obvious guard does not work. Writing

```julia
positive = x > zero(T)
value = positive ? n / x : zero(T)
```

keeps the *primal* safe and is still unsafe under reverse-mode AD: LLVM
speculates `fdiv` out of its guard into a `select`, so the division is
differentiated on the branch whose value is discarded, and at `x == 0` that is
`-n / 0^2 = -Inf`, which the discarded branch's zero cotangent turns into `NaN`.
Substituting a safe denominator inside the guard does not help either -
`positive ? n / one(T) : ...` is folded back to `n / x` because the two are
provably equal wherever the guard holds.

Flooring the denominator unconditionally is what survives that: the division no
longer depends on which branch runs, its value is unchanged wherever
`denominator >= eps(T)`, and `max` contributes a zero derivative below the
floor, so nothing infinite is ever produced. Measured on the maize weather
gradient, this is the difference between a `NaN` gradient and agreement with
finite differences.

Callers keep whatever guard the physics needs for the *value*; this only makes
the arithmetic safe to differentiate.
"""
@inline guarded_quotient(numerator::T, denominator::T) where {T <: AbstractFloat} =
    numerator / max(denominator, eps(T))
