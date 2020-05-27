# Equation references to
#
# - Bethe, H., & Salpeter, E. (1977). Quantum mechanics of one-and
#   two-electron atoms. New York: Plenum Pub. Corp.

@doc raw"""
    non_relativistic_energy(n[, Z=1])

Returns the non-relativistic energy of a hydrogen orbital with
principal quantum number ``n`` (and an effective nuclear charge ``Z``)
according to

```math
\begin{equation}
\tag{BS2.11}
E_n =
-\frac{Z^2}{2n^2}
\end{equation}
```

"""
non_relativistic_energy(n, Z=1) = -Z^2/2n^2

"""
    non_relativistic_energy(o)

Return the non-relativistic energy of the orbital `o`.

# Examples

```jldoctest
julia> non_relativistic_energy(o"1s")
-0.5

julia> non_relativistic_energy(o"2p")
-0.125

julia> non_relativistic_energy(ro"2p")
-0.125

julia> non_relativistic_energy(ro"2p-")
-0.125
```
"""
non_relativistic_energy(o::AbstractOrbital) = non_relativistic_energy(o.n)
non_relativistic_energy(o::SpinOrbital) = non_relativistic_energy(o.orb)

@doc raw"""
    relativistic_energy(n, κ[, Z=1])

Returns the relativistic energy of a hydrogen orbital with principal
quantum number ``n`` and "angular momentum" quantum number ``κ`` (and
an effective nuclear charge ``Z``) according to

```math
\begin{equation}
\tag{BS14.29*}
E_{nκ} =
c^2\left[
-1+
\frac{1}{\sqrt{1 +
\left(
\frac{Zα}{n - \abs{κ} + \sqrt{κ^2-Z^2α^2}}
\right)^2}}
\right]
\end{equation}
```

where the energy has been shifted by ``-c^2`` compared to
``\mathrm{(BS14.29)}``.

"""
function relativistic_energy(n, κ, Z=1)
    r = Z*α/(n-abs(κ)+√(κ^2-Z^2*α^2))
    -c^2*(1-inv(√(1+r^2)))
end

"""
    relativistic_energy(o::RelativisticOrbital)

Return the relativistic energy for the orbital `o`.

# Example

```jldoctest
julia> relativistic_energy(ro"1s")
-0.5000066565961734

julia> relativistic_energy(ro"2p-")
-0.12500208019062972

julia> relativistic_energy(ro"2p")
-0.1250004160284539
```
"""
relativistic_energy(o::RelativisticOrbital) = relativistic_energy(o.n, o.κ)
relativistic_energy(o::SpinOrbital) = relativistic_energy(o.orb)

"""
    orbital_energy(o)

Return the energy of orbital `o`, automatically choosing between
[`non_relativistic_energy`](@ref) and [`relativistic_energy`](@ref),
depending on the type of `o`.

# Examples

```jldoctest
julia> orbital_energy(o"1s")
-0.5

julia> orbital_energy(ro"1s")
-0.5000066565961734

julia> orbital_energy(o"2p")
-0.125

julia> orbital_energy(ro"2p")
-0.1250004160284539
```

"""
orbital_energy(o::Orbital) = non_relativistic_energy(o)
orbital_energy(o::RelativisticOrbital) = relativistic_energy(o)
orbital_energy(o::SpinOrbital) = orbital_energy(o.orb)

ground_state( ::Orbital) = o"1s"
ground_state( ::RelativisticOrbital) = ro"1s"
ground_state(o::SpinOrbital) = ground_state(o.orb)

"""
    atomic_hamiltonian(orbitals[; Efun=orbital_energy, neutral_zero=true])

Return the diagonal matrix whose elements corresponds to the energies
of each of the `orbitals`, by default [`orbital_energy`](@ref) is used
to determine the energy, which dispatches on the orbital type to
return the appropriate energy (non-relativistic or relativistic), and
all energies are shifted such that `1s` has zero energy, if
`neutral_zero==true`.

# Examples

```jldoctest
julia> atomic_hamiltonian(os"1[s] 2[s-p]")
3×3 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.0   ⋅      ⋅
  ⋅   0.375   ⋅
  ⋅    ⋅     0.375

julia> atomic_hamiltonian(ros"1[s] 2[s-p]")
4×4 LinearAlgebra.Diagonal{Float64,Array{Float64,1}}:
 0.0   ⋅         ⋅         ⋅
  ⋅   0.375005   ⋅         ⋅
  ⋅    ⋅        0.375005   ⋅
  ⋅    ⋅         ⋅        0.375006
```
"""
function atomic_hamiltonian(orbitals; Efun::Function=orbital_energy,
                            neutral_zero=true)
    E = zeros(length(orbitals))

    for (i,o) in enumerate(orbitals)
        E[i] = Efun(o)
    end

    if neutral_zero && !isempty(orbitals)
        E .-= Efun(ground_state(first(orbitals)))
    end

    Diagonal(E)
end

export non_relativistic_energy, relativistic_energy, orbital_energy, atomic_hamiltonian
