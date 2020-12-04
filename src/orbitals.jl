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

"""
    radial_moment(k, n, ℓ[, Z = 1])

Compute the radial moment ``⟨rᵏ⟩`` for a hydrogenic orbital ``(n,ℓ)``
in a Coulomb potential of charge ``Z``. The formulas are taken from Table 2⁵ of

- Condon, E. U., & Shortley, G. H. (1951). The theory of atomic
  spectra. Cambridge: Cambridge University Press.

and are available for ``k ∈ -4:4``.
"""
function radial_moment(k::Integer, n::Integer, ℓ::Integer, Z = 1)
    if k == 1
        (3n^2 - ℓ*(ℓ+1))/2Z
    elseif k == 2
        n^2*(5n^2 + 1 - 3ℓ*(ℓ+1))/2Z^2
    elseif k == 3
        n^2*(35n^2*(n^2-1) - 30n^2*(ℓ+2)*(ℓ-1) + 3*(ℓ+2)*(ℓ+1)*ℓ*(ℓ-1))/8Z^3
    elseif k == 4
        n^4*(63n^4 - 35n^2*(2ℓ^2+2ℓ-3) + 5ℓ*(ℓ+1)*(3ℓ^2+3ℓ-10) + 12)/8Z^4
    elseif k == -1
        Z/n^2
    elseif k == -2
        Z^2/(n^3*(ℓ+1//2))
    elseif k == -3
        Z^3/(n^3*(ℓ+1)*(ℓ+1//2)*ℓ)
    elseif k == -4
        Z^4/2*(3n^2 - ℓ*(ℓ+1))/(n^5*(ℓ+3//2)*(ℓ+1)*(ℓ+1//2)*ℓ*(ℓ-1//2))
    elseif k == 0
        one(Z)
    else
        throw(ArgumentError("Radial moment not available for k = $(k)"))
    end
end

export non_relativistic_orbital, radial_moment
