module Hydrogen

using LinearAlgebra
using SparseArrays
using Combinatorics

using AtomicLevels
using AngularMomentumAlgebra
import AngularMomentumAlgebra: powneg1
using AngularMomentumAlgebra.Dipoles
using HypergeometricFunctions
using Polynomials: Polynomial
using SpecialPolynomials: Laguerre
using SpecialFunctions: gamma

include("sparse_builder.jl")

include("constants.jl")
include("energies.jl")
include("dipoles.jl")
include("orbitals.jl")

end
