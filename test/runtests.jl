using StructuredQR
using Test
using BenchmarkTools, GenericLinearAlgebra
using LinearAlgebra
using BlockDiagonals

## Functions useful for tests ##
include("MatricesRebuild.jl")
include("QRutils.jl")

## Tests if the functions work properly ##
include("TestDenseQR.jl") #Successful tests
include("TestBlocDiagQR.jl") #Successful tests
include("TestHCatQR.jl") #Successful tests

## Tests the number of allocations ##
include("AllocsQRDense.jl")
include("AllocsQops.jl")
include("AllocsQRBlocDiag.jl")
include("AllocsQRHcat.jl")

## Benchmarks to show the memory allocations and the elapsed time ##
include("BenchmarkQRDense.jl")
include("BenchmarkQops.jl")
include("BenchmarkQRBlocDiag.jl")
include("BenchmarkQRHcat.jl")