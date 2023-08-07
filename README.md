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
\begin{choices}
m \geq n \\
rank(A) = n
\end{choices}
$$

This matrix $A$ has a unique QR-decomposition if we impose the sign of the diagonal coefficients of $R$. With this factorization, some useful operations with $Q$ can be done (such like the product of $Q$ by a vector or a matrix).

Depending on the type of `A` (and so far its structure as well), a specific QR-factorisation will be executed.

### Dense matrices

For a dense matrix, `qrH!` will be call.

````JULIA
m = 10
n = 8
A = rand(m, n)
qrhat!(A)
````

The QR factorisation is computed inside of `A`, so its previous coefficients have been now replaced by those of $R$ in its upper triangle and the remaining is filled with vectors required to apply operations linked with $Q$.

### Block-Diagonal matrices

For Block-Diagonal matrices defined with the `BlockDiagonals.jl` module already existing in Julia, `QRblocdiag!` will be call.

````JULIA
using BlockDiagonals

A = BlockDiagonal([rand(4,4), rand(3,2)])
qrhat!(A)
````

### Horizontally concatenated matrices

This module also allows the user to apply a QR-factorisation by splitting horizontally the matrix so that the specific QR factorization will consider those two blocks :

````JULIA
A1, A2 = rand(6, 3), rand(6, 2)
A = HcatMatrix(A1, A2)
qrhat!(A)
````

## Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/JuliaSmoothOptimizers/JSOTemplate.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) organization, so questions about any of our packages are welcome.
