"""
    non_relativistic_orbital(n, ℓ[, Z = 1])

Returns a callable object to evaluate the _reduced_ radial orbital ``P_{nℓ}(r)`` of a
hydrogen-like atom in non-relativistic theory, given by

```math
P_{nℓ}(r) =
\\sqrt{
    \\left(\\frac{2Z}{n}\\right)^3
    \\frac{(n-\\ell-1)!}{2n(n+\\ell)!}
}
\\left(\\frac{2Z}{n} \\cdot r\\right)^\\ell
L^{(2\\ell+1)}_{n-\\ell-1}\\left(\\frac{2Z}{n} \\cdot r\\right)
\\exp\\left(-\\frac{Z}{n} \\cdot r\\right)
```

``P_{nℓ}(r)`` is defined through

```math
\\psi_{nlm}(r,\\theta,\\varphi) = \\frac{1}{r} P_{nℓ}(r) Y_{ℓm}(\\theta,\\varphi)
```

where ``\\psi_{nlm}`` are the normalized solutions to the non-relativistic single-particle
Schrödinger equation for a point-like nucleus with charge `Z`, and where ``Y_{ℓm}`` are the
spherical harmonics.
"""
function non_relativistic_orbital(n::Integer, ℓ::Integer, Z = 1)
    n > 0 || throw(ArgumentError("Invalid n = $n, must be n ≥ 1"))
    0 ≤ ℓ < n || throw(ArgumentError("Invalid ℓ = $ℓ for n = $n, must be 0 ≤ ℓ < n"))
    # Normalization factor.
    # Note: n! = Γ(n+1), so (n-ℓ-1)!/(2n(n+ℓ)!) = Γ(n+1)/(2nΓ(n+ℓ+1))
    N = sqrt((2Z/n)^3*gamma(n-ℓ)/(2n*gamma(n+ℓ+1)))
    # Generalized Laguerre polynomial L^{(2ℓ+1)}_{n-ℓ-1}(q), where q = 2Zr/n:
    L = convert(Polynomial, Laguerre{2ℓ+1}([zeros(n-ℓ-1); 1], :q))
    # Single Polynomial object for the polynomial part, multiplied by r to get the
    # reduced polynomial:
    C = let r = Polynomial([0, 1], :r)
        N * (2*Z*r/n)^ℓ * L(2*Z*r/n) * r
    end
    # Multiply with the exponential exp(Zr/n) and return a callable object:
    q = -Z/n # pre-compute the coefficient in the exponential
    r -> C(r) * exp(q*r)
end

export non_relativistic_orbital
