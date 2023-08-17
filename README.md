# StructuredQR package

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/JSOTemplate.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/JSOTemplate.jl/dev)
[![Build Status](https://github.com/JuliaSmoothOptimizers/JSOTemplate.jl/workflows/CI/badge.svg)](https://github.com/JuliaSmoothOptimizers/JSOTemplate.jl/actions)
[![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/JSOTemplate.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/JSOTemplate.jl)
[![Coverage](https://codecov.io/gh/JuliaSmoothOptimizers/JSOTemplate.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSmoothOptimizers/JSOTemplate.jl)

## How to Cite

If you use StructuredQR.jl in your work, please cite using the format given in [CITATION.cff](https://github.com/JuliaSmoothOptimizers/JSOTemplate.jl/blob/main/CITATION.cff).

## Philosophy

Whenever a matrix has a specific well-known structure, taking it into account for factorizations can help to save time and memory allocations. So that, this package has been made to exploit some specific structures of matrices that can be found in some applied problems to compute efficiently a QR factorisation of it.

## Compatibility

Julia 1.8 and up.

## How to Install

````JULIA
pkg> add StructuredQR
pkg> test StructuredQR
````

## How to use

First of all, the QR factorizations presented here have been devised for full-rank overdetermined matrices. It means, for a matrix $A$ of dimension $m \times n$, we got :

$$
\begin{cases}
m \geq n \\
rank(A) = n
\end{cases}
$$

This matrix $A$ has a unique QR-decomposition if we impose the sign of the diagonal coefficients of $R$. With this factorization, some useful operations with $Q$ can be done (such like the product of $Q$ by a vector or a matrix).

Depending on the type of `A` (and so far its structure as well), a specific QR-factorisation will be executed.

### Dense matrices

For a dense matrix (i.e. a matrix without any specific structure), `qrH!` is the function applying the QR-Factorization using Reflections of Householder. It can be used just as follows :

````JULIA
m = 10
n = 8
A = rand(m, n)
qrH!(A)
````

The QR factorisation is computed inside of `A`, so its previous coefficients have been now replaced by those of $R$ in its upper triangle and the remaining is filled with vectors required to apply operations linked with $Q$. Be careful : this function modifies the elements of `A`.

The QR-hat function `qrhat!` can be used instead whenever you have an `AbstractMatrix` in input.

````JULIA
m = 10
n = 8
A = rand(m, n)
qrhat!(A)
````

### Block-Diagonal matrices

For Block-Diagonal matrices defined with the `BlockDiagonals.jl` module already existing in Julia, `QRblocdiag!` can be used to exploit this specific structure.

````JULIA
using BlockDiagonals

A = BlockDiagonal([rand(4,4), rand(3,2)])
QRblocdiag!!(A)
````

Just like for dense matrices, `qrhat!` can be used instead and will apply the suitable method for this structure.

````JULIA
using BlockDiagonals

A = BlockDiagonal([rand(4,4), rand(3,2)])
qrhat!(A)
````

### Horizontally concatenated matrices

To treat the case of horizontally concatenated matrices, this module uses the structure of `BlockArray` from the module `BlockArrays.jl` but only the very particular case of a 1x2-Block-Matrix.

````JULIA
using BlockArrays

A = BlockArray(rand(6,5), [6], [2,3])
QRhcat!(A)
````

Like the two previous cases, `qrhat!` can be used to apply the suitable method for horizontally concatenated matrices whenever the function has in input a Block-Array.

````JULIA
using BlockArrays

A = BlockArray(rand(6,5), [6], [2,3])
qrhat!(A)
````

## Specific features for this package

One of the relevant specific features of this package is that all of the functions available do not allocate additionnal memory. This feature is a real advantage to avoid some `Out of Memory errors` when you handle very large matrices.

A second convenient feature is the hat function devised in this package is here to apply the most suitable method for the structure of the input matrix. But if the user prefers to use directly the most appropriate function for their case, they can simply call the desired QR-factorization function.

Another feature is the `QOperations.jl` in which you have several linear operations using $Q$, the unitary matrix from the QR decomposition. Just like all of the other functions, these functions do not allocate additionnal memory. However, to use these functions, you absolutely need to compute first a QR factorisation just to get the necessary information about Q. But afterwards, you'll can carry out some basic linear operations ($Qb$, $Q^Tb$, $QB$) etc... 

In `QOperations.jl`, there is as well operations to solve some linear systems. The function `qrsolve!` solves the system $Ax = b$ for an overdetermined full-rank matrix $A$ via the QR-Factorization of $A$.

````JULIA
A = rand(10, 8)
b = rand(10)
qrsolve!(A, b)
````
Besides, `rsolve!` solves the system $Rx = y$ where $R$ is an upper triangular matrix.

````JULIA
A = rand(10, 8)
b = rand(10)
qrhat!(A)
rdiv!(A, b)
````

## Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/JSOTemplate.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) organization, so questions about any of our packages are welcome.
