module Hydrogen

using LinearAlgebra
using SparseArrays

using AtomicLevels
using AngularMomentumAlgebra
import AngularMomentumAlgebra: powneg1
using AngularMomentumAlgebra.Dipoles
using HypergeometricFunctions

include("sparse_builder.jl")

include("constants.jl")
include("energies.jl")

end
