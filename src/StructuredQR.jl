module StructuredQR

using LinearAlgebra, Printf, SparseArrays, BlockDiagonals, BlockArrays

include("QOperations.jl")
include("QRBlocDiag.jl")
include("QRDense.jl")
include("QRhcat.jl")
include("qrhat.jl")

end
