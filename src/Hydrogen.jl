module Hydrogen

using LinearAlgebra
using SparseArrays

using AtomicLevels
using AngularMomentumAlgebra
import AngularMomentumAlgebra: powneg1
using AngularMomentumAlgebra.Dipoles
using HypergeometricFunctions

include("sparse_builder.jl")
include("factorial_ratio.jl")

include("constants.jl")
include("energies.jl")
include("dipoles.jl")

end
