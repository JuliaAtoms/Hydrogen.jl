#=
Equation numbers refer to

- Burgess, A. (1965). Tables of hydrogenic photoionization
  cross-sections and recombination coefficients. Monthly Notices of
  the Royal Astronomical Society, 69, 1.

=#

function ℓℓ′2index(ℓ,ℓ′)
    if ℓ == 0
        ℓ′ == 1 || throw(ArgumentError("Invalid ℓ′ = $(ℓ′) for ℓ = 0"))
        return 1
    end
    abs(ℓ′-ℓ) == 1 || throw(ArgumentError("Invalid ℓ′ = $(ℓ′) for ℓ = $(ℓ)"))
    2ℓ + (ℓ′>ℓ)
end

struct BurgessG{T}
    n::Int
    g::Matrix{T}
end

function BurgessG(::Type{T}, n::Integer, nκ) where T
    ℓs = 0:n-1
    num_channels = 2length(ℓs)-1
    g = zeros(T, num_channels, nκ)
    BurgessG(n, g)
end

Base.getindex(g::BurgessG, ℓ, ℓ′, iκ) =
    g.g[ℓℓ′2index(ℓ,ℓ′), iκ]

Base.setindex!(g::BurgessG, v, ℓ, ℓ′, iκ) =
    setindex!(g.g, v, ℓℓ′2index(ℓ,ℓ′), iκ)

function recur!(g::BurgessG{T}, κ) where T
    n = g.n
    n == 1 && return

    n² = T(n^2)
    κ² = T.(κ.^2)

    for (j,κ²) in enumerate(κ²)
        a = one(T) + n²*κ²

        # Equation (32)
        g[n-2,n-1,j] =
            √((2n-1)*a) *
            g[n-1,n,j]/T(2)

        # Equation (33)
        g[n-1,n-2,j] =
            √(a/(1 + (n-1)^2*κ²)) *
            g[n-1,n,j]/T(2n)

        if n > 2
            # Equation (34)
            g[n-2,n-3,j] =
                (4+(n-1)*a)/T(2n) *
                √(T(2n-1)/(1 + (n-2)^2*κ²)) *
                g[n-1,n-2,j]
        end
    end

    n == 2 && return
    
    for ℓ = n-1:-1:2
        ℓ² = T(ℓ^2)
        c₁ = n² - ℓ²
        for (j,κ²) in enumerate(κ²)
            c₂ = T(1 + ℓ²*κ²)
            c₃ = T(1 + n²*κ²)

            # Equation (28)
            g[ℓ-2,ℓ-1,j] =
                ((4*c₁ + ℓ*(2ℓ-1)*c₃) *
                g[ℓ-1,ℓ,j] -
                2n*√(c₁*(1+(ℓ+1)^2*κ²)) *
                g[ℓ,ℓ+1,j])/(2n*√((n²-(ℓ-1)^2)*c₂))

            if ℓ+1<n
                # Equation (29)
                g[ℓ-1,ℓ-2,j] =
                    ((4*c₁ + ℓ*(2ℓ+1)*c₃) *
                    g[ℓ,ℓ-1,j] -
                    2n*√((n² - (ℓ+1)^2)*c₂) *
                    g[ℓ+1,ℓ,j])/(2n*√(c₁*(1+(ℓ-1)^2*κ²)))
            end
        end
    end

end

@doc raw"""
    burgess_g(n, κ)

Compute the integrals
```math
\begin{equation}
\tag{Burgess22}
g(n, \ell, \kappa, \ell') =
\frac{1}{n^2}
\int_0^\infty\diff{\rho}
\mathscr{P}_{n\ell}(\rho)
\rho
\mathscr{F}_{\kappa\ell}(\rho)
\end{equation}
```
as given by [Burgess1965](@citet), for ``\ell\in\{0..n-1\}`` and
``\ell'=\ell\pm1``.
"""
function burgess_g(n::Integer, κ::Union{Real,AbstractVector{<:Real}})
    T = eltype(κ)
    κ′ = vcat(zero(T), κ)

    g = BurgessG(T, n, length(κ′))

    # Equation (30)
    g[n-1,n,1] = √(π/(2*factorial(2n-1)))*4*(4n)^n*exp(-2n)

    for (i,κ) ∈ enumerate(κ′)
        i == 1 && continue
        if iszero(κ)
            copyto!(view(g.g, :, i), view(g.g, :, 1))
            continue
        end
        κ² = κ^2

        # Equation (31)
        a = √(prod(s -> 1+s^2*κ², 1:n, init=one(T))/(1-exp(-2π/κ)))
        b = exp(2n-2inv(κ)*atan(n*κ))/(1+n^2*κ²)^(n+2)
        g[n-1,n,i] = a*b*g[n-1,n,1]
    end

    recur!(g, κ′)

    g.g[:,2:end]
end
