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
\frac{Zα}{n - |κ| + \sqrt{κ^2-Z^2α^2}}
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
