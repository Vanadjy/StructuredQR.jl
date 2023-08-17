using StructuredQR
using Test
using TestSetExtensions
using GenericLinearAlgebra
using LinearAlgebra
using BlockDiagonals
using BlockArrays

## Functions useful for tests ##
include("MatricesRebuild.jl")
include("QRutils.jl")

## Tests if the functions work properly ##
include("TestDenseQR.jl") 
include("TestBlocDiagQR.jl") 
include("TestHCatQR.jl")
include("TestQRsolve.jl")

## Tests the number of allocations ##
println("-------------------------------")
println("Testing allocations")
println("-------------------------------")

@testset "Allocations of functions" begin
include("AllocsQRDense.jl")
include("AllocsQops.jl")
include("AllocsQRBlocDiag.jl")
include("AllocsQRHcat.jl")
include("AllocsQRhat.jl")
include("AllocsQRsolve.jl")
end