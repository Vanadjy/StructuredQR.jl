export qprod!, qtprod!, qprod, qtprod, qmul!, qtmul!, qmul, qtmul, get_r, rdiv!, rdiv, qrsolve!, q_lop

"""
qtprod!(A::AbstractMatrix,x::AbstractVector)

Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

Computes : Q*x

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `x`: a vector of size m containing the result of the operation Q*x;
"""

function qtprod!(A::AbstractMatrix,x::AbstractVector)
    m, n = size(A)
    j = 1
    while (j <= n) & (j < m)
        uj = view(A,j+1:m,j)
        δj = uj'uj + 1
        xⱼ = x[j]
        β = 0
        for i = 1:m-j
            β += uj[i]*x[i+j]
        end
        x[j] -= 2*(xⱼ + β)/δj
        for l = j+1:m
            x[l] -= 2*(xⱼ + β)*A[l,j]/δj
        end
        j+=1
    end
    x
end

function qtprod!(A::BlockDiagonal{T, Matrix{T}}, x::AbstractVector) where T
    j = 1
    N = nblocks(A)
    indexi = 0
    while j ≤ N
        m = size(BlockDiagonals.blocks(A)[j], 1)
        @views qtprod!(BlockDiagonals.blocks(A)[j], x[(indexi + 1):(indexi + m)])
        indexi += m
        j += 1
    end
end

"""
qprod!(A::AbstractMatrix,x::AbstractVector)

Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

Computes : Qx

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `x`: a vector of size m containing the result of the operation Qx;
"""
function qprod!(A::AbstractMatrix,x::AbstractVector)
    m, n = size(A)
    k = 1
    if m==n #condition to treat the square case
        k +=1
    end
    while (k <= n) & (k <= m)
        j = n-k+1
        uj = view(A,j+1:m,j)
        δj = uj'uj + 1
        xⱼ = x[j]
        β = 0
        for i = 1:m-j
            β += uj[i]*x[i+j]
        end
        x[j] -= 2*(xⱼ + β)/δj
        for l = j+1:m
            x[l] -= 2*(xⱼ + β)*A[l,j]/δj
        end
        k+=1
    end
    x
end

function qprod!(A::BlockDiagonal{T, Matrix{T}}, x::AbstractVector) where T
    j = 1
    N = nblocks(A)
    indexi = 0
    while j ≤ N
        m = size(BlockDiagonals.blocks(A)[j], 1)
        @views qprod!(BlockDiagonals.blocks(A)[j], x[(indexi + 1):(indexi + m)])
        indexi += m
        j += 1
    end
end

"""
qprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector)

Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

Computes : Qx

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `y`: a newly created vector of size m containing the result of the operation Qx;
"""

function qprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector) where T
    y = copy(x)
    qprod!(A, y)
    y
end

"""
qtprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector)

Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

Computes : Q*x

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `y`: a newly created vector of size m containing the result of the operation Q*x;
"""

function qtprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector) where T
    y = copy(x)
    qtprod!(A, y)
    y
end

"""
qmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QB

Computes : QB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `B`: a matrix of dimensions m × r containing the result of the operation QB;
"""

function qmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    m, n = size(B)
    for j = 1:n
        #we apply here qprod! on each column of B : Calculates Qbⱼ
        qprod!(A, view(B, 1:m, j))
    end
    B
end

"""
qmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within Y by replacing its values by those of QB

Computes : QB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `Y`: a newmy created matrix of dimensions m × r containing the result of the operation QB;
"""

function qmul(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    Y = copy(B)
    qmul!(A, Y)
    Y
end

"""
qtmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QᵀB

Computes : QᵀB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `B`: a matrix of dimensions m × r containing the result of the operation QᵀB;
"""

function qtmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    m, n = size(B)
    for j = 1:n
        #we apply here qtprod! on each column of B : Calculates Qᵀbⱼ
        qtprod!(A, view(B, 1:m, j))
    end
    B
end

"""
qtmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within Y by replacing its values by those of QᵀB

Computes : QᵀB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `Y`: a newly created matrix of dimensions m × r containing the result of the operation QᵀB;
"""

function qtmul(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    Y = copy(B)
    qtmul!(A, Y)
    Y
end

"""
get_r(A::Union{AbstractMatrix, AbstractSparseMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}})

Returns the R factor from the QR factorization of A

#### Input arguments
* `A` : a full-rank overdetermined matrix of dimension (m × n) that has been QR-factorized

#### Output arguments
* `R` : the R factor of dimension (m × n) (uppper triangular matrix) of the QR-factorization of A
"""

## Adding methods for Linear Operator ##

"""
qtprod!(res::AbstractVector, A::AbstractMatrix,x::AbstractVector)

Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

Computes : Q*x

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 
* `res`: a vector similar to x

#### Output arguments

* `res`: a vector of size m containing the result of the operation Q*x;
"""

function qtprod!(res::AbstractVector, A::AbstractMatrix, x::AbstractVector)
    m, n = size(A)
    res .= x
    j = 1
    while (j <= n) & (j < m)
        uj = view(A,j+1:m,j)
        δj = uj'uj + 1
        xⱼ = res[j]
        β = 0
        for i = 1:m-j
            β += uj[i]*res[i+j]
        end
        res[j] -= 2*(xⱼ + β)/δj
        for l = j+1:m
            res[l] -= 2*(xⱼ + β)*A[l,j]/δj
        end
        j+=1
    end
    res
end

function qtprod!(res::AbstractVector, A::BlockDiagonal{T, Matrix{T}}, x::AbstractVector) where T
    j = 1
    N = nblocks(A)
    indexi = 0
    while j ≤ N
        m = size(BlockDiagonals.blocks(A)[j], 1)
        @views qtprod!(res[(indexi + 1):(indexi + m)], BlockDiagonals.blocks(A)[j], x[(indexi + 1):(indexi + m)])
        indexi += m
        j += 1
    end
    res
end

"""
qprod!(res::AbstractVector, A::AbstractMatrix,x::AbstractVector)

Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

Computes : Qx

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 
* `res`: a vector similar to x

#### Output arguments

* `res`: a vector of size m containing the result of the operation Qx;
"""
function qprod!(res::AbstractVector, A::AbstractMatrix,x::AbstractVector)
    m, n = size(A)
    res .= x
    k = 1
    if m==n #condition to treat the square case
        k +=1
    end
    while (k <= n) & (k <= m)
        j = n-k+1
        uj = view(A,j+1:m,j)
        δj = uj'uj + 1
        xⱼ = res[j]
        β = 0
        for i = 1:m-j
            β += uj[i]*res[i+j]
        end
        res[j] -= 2*(xⱼ + β)/δj
        for l = j+1:m
            res[l] -= 2*(xⱼ + β)*A[l,j]/δj
        end
        k+=1
    end
    res
end

function qprod!(res::AbstractVector, A::BlockDiagonal{T, Matrix{T}}, x::AbstractVector) where T
    j = 1
    N = nblocks(A)
    indexi = 0
    while j ≤ N
        m = size(BlockDiagonals.blocks(A)[j], 1)
        @views qprod!(res[(indexi + 1):(indexi + m)], BlockDiagonals.blocks(A)[j], x[(indexi + 1):(indexi + m)])
        indexi += m
        j += 1
    end
end

"""
qprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector)

Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

Computes : Qx

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `res`: a newly created vector of size m containing the result of the operation Qx;
"""

function qprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector) where T
    res = copy(x)
    qprod!(res, A , x)
    res
end

"""
qtprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector)

Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

Computes : Q*x

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - x is a vector

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `x`: a vector of size m 

#### Output arguments

* `res`: a newly created vector of size m containing the result of the operation Q*x;
"""

function qtprod(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}},x::AbstractVector) where T
    res = copy(x)
    qtprod!(res, A, x)
    res
end

"""
qmul!(res::AbstractMatrix , A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T

Calculates the multiplication of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QB

Computes : QB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r
* `res`: a matrix similar to B

#### Output arguments

* `res`: a matrix of dimensions m × r containing the result of the operation QB;
"""

function qmul!(res::AbstractMatrix , A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    m, n = size(B)
    for j = 1:n
        #we apply here qprod! on each column of B : Calculates Qbⱼ
        qprod!(view(res, 1:m, j), A, view(B, 1:m, j))
    end
    B
end

"""
qmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within Y by replacing its values by those of QB

Computes : QB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `res`: a newmy created matrix of dimensions m × r containing the result of the operation QB;
"""

function qmul(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    res = copy(B)
    qmul!(res, A, B)
    res
end

"""
qtmul!(res::AbstractMatrix ,A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T

Calculates the multiplication of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QᵀB

Computes : QᵀB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r
* `res`: a matrix similar to B

#### Output arguments

* `res`: a matrix of dimensions m × r containing the result of the operation QᵀB;
"""

function qtmul!(res::AbstractMatrix, A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    m, n = size(B)
    for j = 1:n
        #we apply here qtprod! on each column of B : Calculates Qᵀbⱼ
        @views qtprod!(res[1:m, j], A, B[1:m, j])
    end
    res
end

"""
qtmul!(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix)

Calculates the multiplication of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within Y by replacing its values by those of QᵀB

Computes : QᵀB

Where :
    - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
    - B is a matrix

#### Input arguments

* `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
* `B`: a matrix of dimensions m × r

#### Output arguments

* `res`: a newly created matrix of dimensions m × r containing the result of the operation QᵀB;
"""

function qtmul(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}, B::AbstractMatrix) where T
    res = copy(B)
    qtmul!(res, A, B)
    res
end

function q_lop(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}}) where T
    nrow = size(A)[1]
    #multiplication of x by Q which is squared
    op_qprod = LinearOperator(AbstractFloat, nrow, nrow, false, false,
                            (res, x) -> qprod!(res, A, x),
                            nothing,
                            (res, y) -> qtprod!(res, A, y))
    return op_qprod
end

"""
get_r(A::Union{AbstractMatrix, AbstractSparseMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}})

Returns the R factor from the QR factorization of A

#### Input arguments
* `A` : a full-rank overdetermined matrix of dimension (m × n) that has been QR-factorized

#### Output arguments
* `R` : the R factor of dimension (m × n) (uppper triangular matrix) of the QR-factorization of A
"""
function get_r(A::Union{AbstractMatrix, AbstractSparseMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}}) where T
    return triu(A)
end

"""
rdiv!(res::AbstractVector, A::Union{AbstractMatrix, AbstractBlockMatrix{T}}, b::AbstractVector) where T

    Computes the solution of the system Rx = b where R is upper triangular by replacing the values of the rhs b

    #### Input arguments
    * `A` : a full-rank overdetermined matrix of dimension (m × n) that has been QR-factorized
    * `b` : a vector of size m
    * `res`: a vector similar to b

    #### Output arguments
    `res` : the vector containing the solution of the system Rx = b     
"""
function rdiv!(res::AbstractVector, A::Union{AbstractMatrix, AbstractBlockMatrix{T}}, b::AbstractVector) where T
    n = size(A, 2)
    res .= b
    for j = n:-1:1
        @views res[j] = (res[j] - A[j, (j + 1):n]'res[(j + 1):n])/A[j,j]
    end
    res
end

"""
    rdiv!(res::AbstractVector, A::BlockDiagonal{T, Matrix{T}}, b::AbstractVector) where T

    Computes the solution of the system Ax = b where A is BlockDiagonal upper triangular by replacing the values of the rhs b. For A with r blocks, r systems are solved using a dense method of rdiv! on each. 

    #### Input arguments
    * `A` : a full-rank overdetermined Block-Diagonal matrix of dimension (m × n) that has been QR-factorized
    * `b` : a vector of size m
    * `res`: a vector similar to b

    #### Output arguments
    `res` : the vector containing the solution of the system Rx = b      
"""
function rdiv!(res::AbstractVector, A::BlockDiagonal{T, Matrix{T}}, b::AbstractVector) where T
    j = 1
    N = nblocks(A)
    indexi = 0
    while j ≤ N
        m = size(BlockDiagonals.blocks(A)[j], 1)
        @views rdiv!(res[(indexi + 1):(indexi + m)], BlockDiagonals.blocks(A)[j], b[(indexi + 1):(indexi + m)])
        indexi += m
        j += 1
    end
    res
end

function rdiv(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}}, b::AbstractVector) where T
    res = copy(b)
    rdiv!(res, A, b)
    res
end

"""
qrsolve!(res::AbstractVector, A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}}, b::AbstractVector) where T

Solves the system Ax = b via the QR factorization of A, so that it solves Rx = Qᵀb

#### Input arguments
`A` : a full-rank overdetermined matrix of dimension (m × n)

#### Output arguments
`b` : the rhs containing now the solution of the system Ax = b

#### References

* TREFETHEN, Lloyd N. et BAU, David. Numerical linear algebra. Siam, 2022, (pp. 54).
"""
function qrsolve!(res::AbstractVector, A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}}, b::AbstractVector) where T
    #Computes the QR factorization
    qrhat!(A)
    #Transforms the rhs as Qᵀb
    qtprod!(res, A, b)
    #Solves the system Rx = Qᵀb
    rdiv!(res, A, res)
end

function qrsolve(A::Union{AbstractMatrix, BlockDiagonal{T, Matrix{T}}, AbstractBlockMatrix{T}}, b::AbstractVector) where T
    res = copy(b)
    qrsolve!(res, A, b)
    res
end

m, n = 10, 8
A = rand(m, n)
B = deepcopy(A)
qtmul!(B, A, B)
Res = similar(B)
a = @allocated begin
    qtmul!(B, A, B)
end; show(a)