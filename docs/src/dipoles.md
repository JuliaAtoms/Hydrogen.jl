# Dipoles

```@docs
Hydrogen.dipole_moment
Hydrogen.radial_dipole_moment
Hydrogen.dipole_matrix
```

## Photoionization dipole matrix elements

We implement the recurrence relations for computing the radial dipole
integrals between a bound orbital and a continuum orbital, as detailed
by [Burgess1965](@citet). NB that his derivation is in _Rydberg atomic
units_, whereas we are typically interested in obtaining results in
_Hartree atomic units_. To reduce the possibilities for errors, the
implementation is in Rydberg atomic units, following the notation of
[Burgess1965](@citet), and an interface routine is provided that
performs unit conversions.

```@docs
Hydrogen.burgess_g
```
