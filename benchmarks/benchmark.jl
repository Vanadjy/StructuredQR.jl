## Benchmarks to show the memory allocations and the elapsed time ##
using StructuredQR
using BenchmarkTools
using GenericLinearAlgebra
using LinearAlgebra
using BlockDiagonals
using BlockArrays

include("CreateDiagBlock.jl")
include("BenchmarkQRDense.jl")
include("BenchmarkQops.jl")
include("BenchmarkQRBlocDiag.jl")
include("BenchmarkQRHcat.jl")
include("BenchmarkQRhat.jl")
include("BenchmarkQRsolve.jl")