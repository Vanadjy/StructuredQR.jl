module StructuredQR

using LinearAlgebra, Printf, SparseArrays, BlockDiagonals, BlockArrays, LinearOperators

include("QRBlocDiag.jl")
include("QRDense.jl")
include("QRhcat.jl")
include("qrhat.jl")
include("QOperations.jl")

end
