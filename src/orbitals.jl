"""
    non_relativistic_orbital(n, ℓ[, Z = 1])

Returns a callable object to evaluate the _reduced_ radial orbital ``P_{nℓ}(r)`` of a
hydrogen-like atom in non-relativistic theory.

``P_{nℓ}(r)`` is defined through

```math
\\psi_{nlm}(r,\\theta,\\varphi) = \\frac{1}{r} P_{nℓ}(r) Y_{ℓm}(\\theta,\\varphi)
```

where ``\\psi_{nlm}`` are the solutions to the non-relativistic single-particle Schrödinger
equation for a point-like nucleus with charge `Z`, and ``Y_{ℓm}`` are the spherical
harmonics.
"""
function non_relativistic_orbital(n::Integer, ℓ::Integer, Z = 1)
    n > 0 || throw(ArgumentError("Invalid n = $n, must be n ≥ 1"))
    0 ≤ ℓ < n || throw(ArgumentError("Invalid ℓ = $ℓ for n = $n, must be 0 ≤ ℓ < n"))
    # Construct a Polynomial object for the polynomial part of P:
    r = Polynomial([0, 1], :r)
    # Normalization factor.
    # Note: n! = Γ(n + 1), so (n - ℓ - 1)!/(2n(n+ℓ)) = Γ(n + 1)/(2nΓ(n+ℓ+1))
    N = sqrt((2Z/n)^3*gamma(n-ℓ)/(2n*gamma(n+ℓ+1)))
    # Generalized Laguerre polynomial L^{(2ℓ+1)}_{n - ℓ - 1}(q), where q = 2Zr/n:
    Lq = convert(Polynomial, Laguerre{2ℓ+1}([zeros(n - ℓ - 1); 1], :q))
    # Single Polynomial object for the polynomial part, multiplied by r to get the
    # reduced polynomial:
    Cr = N * (2*Z*r/n)^ℓ * Lq(2*Z*r/n) * r
    # Multiply with the exponential exp(Zr/n) and return a callable object:
    r -> Cr(r) * exp(-Z*r/n)
end

export non_relativistic_orbital
