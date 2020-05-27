@doc raw"""
    radial_dipole_moment(a, b)

Compute the radial dipole moment between the orbitals `a` and `b` according to

```math
\begin{align}
\tag{BS63.2}
R_{n \ell}^{n',\ell-1} &=
ABCD,
&
n \neq n'\\
\tag{BS63.5}
R_{n,\ell-1}^{n\ell} &= 
R_{n\ell}^{n,\ell-1} =
\frac{3}{2}
n\sqrt{n^2-\ell^2},
& n = n',
\end{align}
```
where
```math
\begin{aligned}
A &= \frac{(-)^{n'-\ell}}{4(2\ell-1)!}, \\
B &= \sqrt{\frac{(n+\ell)!(n'+\ell-1)!}{(n-\ell-1)!(n'-\ell)!}}, \\
C &= \frac{(4nn')^{\ell+1}(n-n')^{n+n'-2\ell-2}}{(n+n')^{n+n'}}, \\
D &= {_2}F_1(-n_r, -n_r', 2\ell, x) -
\left(\frac{n-n'}{n+n'}\right)^2
{_2}F_1(-n_r-2, -n_r', 2\ell, x), \\
n_r &= n-\ell-1, \quad
n_r' = n' - \ell, \quad
x = -\frac{4nn'}{(n-n')^2},
\end{aligned}
```
and ``{_2}F_1(a,b,c;x)`` is a hypergeometric function that can be
computed used
[HyperGeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl).

"""
function radial_dipole_moment(a::Orbital, b::Orbital)
    a == b && return 0.0

    n,ℓ = a.n,a.ℓ
    n′,ℓ′ = b.n,b.ℓ

    abs(ℓ-ℓ′) ≠ 1 && return 0.0
    ℓ < ℓ′ && return radial_dipole_moment(b, a)
    
    # Eq. (BS63.5)
    n == n′ && return 3/2*n*√(n^2-ℓ^2)

    A = powneg1(n′-ℓ)/(4*factorial(2ℓ-1))
    B = √(factorial_ratio(n+ℓ,n-ℓ-1)*factorial_ratio(n′+ℓ-1,n′-ℓ))
    C = (4*n*n′)^(ℓ+1)*(n-n′)^(n+n′-2ℓ-2)/((n+n′)^(n+n′))

    nᵣ = n-ℓ-1
    nᵣ′ = n′-ℓ
    x = -(4n*n′)/((n-n′)^2)
    D = _₂F₁(-nᵣ, -nᵣ′, 2ℓ, x) - ((n-n′)/(n+n′))^2*_₂F₁(-nᵣ-2, -nᵣ′, 2ℓ, x)

    A*B*C*D
end

export radial_dipole_moment
