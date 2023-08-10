using StructuredQR
using Test
using BenchmarkTools, GenericLinearAlgebra
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

## Tests the number of allocations ##
include("AllocsQRDense.jl")
include("AllocsQops.jl")
include("AllocsQRBlocDiag.jl")
include("AllocsQRHcat.jl")
include("AllocsQRhat.jl")

## Benchmarks to show the memory allocations and the elapsed time ##
include("BenchmarkQRDense.jl")
include("BenchmarkQops.jl")
include("BenchmarkQRBlocDiag.jl")
include("BenchmarkQRHcat.jl")
include("BenchmarkQRhat.jl")