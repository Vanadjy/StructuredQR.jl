module StructuredQR

using LinearAlgebra, Printf, SparseArrays, BlockDiagonals

include("Structures/SpecialStructures.jl")
include("QOperations.jl")
include("QRBlocDiag.jl")
include("QRDense.jl")
include("QRhcat.jl")
include("qrhat.jl")

end
