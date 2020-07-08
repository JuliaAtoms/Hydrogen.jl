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
    B = √(factorial(n+ℓ,n-ℓ-1)*factorial(n′+ℓ-1,n′-ℓ))
    C = (4*n*n′)^(ℓ+1)*(n-n′)^(n+n′-2ℓ-2)/((n+n′)^(n+n′))

    nᵣ = n-ℓ-1
    nᵣ′ = n′-ℓ
    x = -(4n*n′)/((n-n′)^2)
    D = _₂F₁(-nᵣ, -nᵣ′, 2ℓ, x) - ((n-n′)/(n+n′))^2*_₂F₁(-nᵣ-2, -nᵣ′, 2ℓ, x)

    A*B*C*D
end

radial_dipole_moment(a::SpinOrbital, b::SpinOrbital) =
    radial_dipole_moment(a.orb, b.orb)

"""
    dipole_moment(a, b[, component=:z])

Compute the dipole moment between orbitals `a` and `b` along the
Cartesian direction indicated by `component`. The radial integral is
computed using [`radial_dipole_moment`](@ref) and the angular integral
using
[AngularMomentumAlgebra.jl](https://github.com/JuliaAtoms/AngularMomentumAlgebra.jl).

```jldoctest
julia> a,b,c,d = SpinOrbital(o"1s", 0, half(1)), SpinOrbital(o"2s", 0, half(1)), SpinOrbital(o"2p", 0, half(1)), SpinOrbital(o"2p", 1, half(1))
(1s₀α, 2s₀α, 2p₀α, 2p₁α)

julia> dipole_moment(a, c)
0.7449355390278029

julia> dipole_moment(b, c)
2.999999999999999

julia> dipole_moment(a, d)
0.0

julia> dipole_moment(a, d, :x)
-0.5267489711934153
```

"""
function dipole_moment(a::SpinOrbital, b::SpinOrbital, component=:z)
    d = AngularMomentumAlgebra.Dipoles.𝐫[(x = 1, y = 2, z = 3)[component]]
    v = dot(a,d,b)
    iszero(v) && return 0.0
    r = radial_dipole_moment(a, b)
    r*v.operators[1][2]
end

"""
    dipole_matrix(orbitals[, component=:z])

Return the sparse matrix whose elements correspond to the dipole
moments between all of the `orbitals` as measured along the Cartesian
direction indicated by `component`.

# Examples

```jldoctest
julia> dipole_matrix(sos"1[s] 2[s-p]")
10×10 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:
  [7, 1]  =  0.744936
  [8, 2]  =  0.744936
  [7, 3]  =  3.0
  [8, 4]  =  3.0
  [1, 7]  =  0.744936
  [3, 7]  =  3.0
  [2, 8]  =  0.744936
  [4, 8]  =  3.0

julia> dipole_matrix(sos"1[s] 2[s-p]", :y)
10×10 SparseArrays.SparseMatrixCSC{Complex{Float64},Int64} with 16 stored entries:
  [5 ,  1]  =  0.0+0.526749im
  [9 ,  1]  =  0.0+0.526749im
  [6 ,  2]  =  0.0+0.526749im
  [10,  2]  =  0.0+0.526749im
  [5 ,  3]  =  0.0+2.12132im
  [9 ,  3]  =  0.0+2.12132im
  [6 ,  4]  =  0.0+2.12132im
  [10,  4]  =  0.0+2.12132im
  [1 ,  5]  =  -0.0-0.526749im
  [3 ,  5]  =  -0.0-2.12132im
  [2 ,  6]  =  -0.0-0.526749im
  [4 ,  6]  =  -0.0-2.12132im
  [1 ,  9]  =  -0.0-0.526749im
  [3 ,  9]  =  -0.0-2.12132im
  [2 , 10]  =  -0.0-0.526749im
  [4 , 10]  =  -0.0-2.12132im
```
"""
function dipole_matrix(orbitals::AbstractVector{<:SpinOrbital}, component=:z)
    d = AngularMomentumAlgebra.Dipoles.𝐫[(x = 1, y = 2, z = 3)[component]]

    T = component == :y ? ComplexF64 : Float64
    norb = length(orbitals)

    sparse_builder(T, norb, norb) do set!
        for (i,a) in enumerate(orbitals)
            for (j,b) in enumerate(orbitals)
                v = dot(a,d,b)
                iszero(v) && continue
                r = radial_dipole_moment(a.orb, b.orb)
                set!(i, j, r*v.operators[1][2])
            end
        end
    end
end

export radial_dipole_moment, dipole_moment, dipole_matrix
