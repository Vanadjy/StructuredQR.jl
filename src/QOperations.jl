export qprod!, qtprod!, qprod, qtprod, qmul!, qtmul!

function qtprod!(A::AbstractMatrix,x::AbstractVector)
    """
    qtprod!(A::AbstractMatrix,x::AbstractVector)

    Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

    Comutes : Q*x

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - x is a vector

    #### Input arguments

    * `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `x`: a vector of size m 
    
    #### Output arguments
    
    * `x`: a vector of size m containing the result of the operation Q*x;
    """
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

function qprod!(A::AbstractMatrix,x::AbstractVector)
    """
    qprod!(A::AbstractMatrix,x::AbstractVector)

    Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within x by replacing its values by those of Qᵀx

    Comutes : Qx

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - x is a vector

    #### Input arguments

    * `A`: an overtermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `x`: a vector of size m 
    
    #### Output arguments
    
    * `x`: a vector of size m containing the result of the operation Qx;
    """
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

function qprod(A::AbstractMatrix,x::AbstractVector)
    """
    qprod(A::AbstractMatrix,x::AbstractVector)

    Calculates the product of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

    Comutes : Qx

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - x is a vector

    #### Input arguments

    * `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `x`: a vector of size m 
    
    #### Output arguments
    
    * `y`: a newly created vector of size m containing the result of the operation Qx;
    """
    y = similar(x, size(A, 1))
    qprod!(A, y)
    y
end

function qtprod(A::AbstractMatrix,x::AbstractVector)
    """
    qtprod(A::AbstractMatrix,x::AbstractVector)

    Calculates the product of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by a vector x and stores the result within y, a new vector with the same size as x

    Comutes : Q*x

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - x is a vector

    #### Input arguments

    * `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `x`: a vector of size m 
    
    #### Output arguments
    
    * `y`: a newly created vector of size m containing the result of the operation Q*x;
    """
    y = similar(x, size(A, 1))
    qtprod!(A, y)
    y
end

function qmul!(A::AbstractMatrix, B::AbstractMatrix)
    """
    qmul!(A::AbstractMatrix, B::AbstractMatrix)

    Calculates the multiplication of Q, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QB

    Comutes : QB

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - B is a matrix

    #### Input arguments

    * `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `B`: a matrix of dimensions m × r
    
    #### Output arguments
    
    * `B`: a matrix of dimensions m × r containing the result of the operation QB;
    """
    m, n = size(B)
    for j = 1:n
        #on effectue ici à chaque fois le produit compact Q*bⱼ
        qprod!(A, view(B, 1:m, j)) #produit vecteur colonne par vecteur colonne de B
    end
    B
end

function qtmul!(A::AbstractMatrix, B::AbstractMatrix)
    """
    qtmul!(A::AbstractMatrix, B::AbstractMatrix)

    Calculates the multiplication of Qᵀ, the unitary matrix from the QR decomposition of A (an overdetermined full-rank matrix), by an other matrix B and stores the result within it by replacing its values by those of QᵀB

    Comutes : QᵀB

    Where :
        - Q is the unitary matrix from the QR factorization of A (an overdetermined full-rank matrix)
        - B is a matrix

    #### Input arguments

    * `A`: an overdetermined full rank matrix of dimension m × n, m ≥ n, containing the coefficients of Q and R;
    * `B`: a matrix of dimensions m × r
    
    #### Output arguments
    
    * `B`: a matrix of dimensions m × r containing the result of the operation QᵀB;
    """
    m, n = size(B)
    for j = 1:n
        #on effectue ici à chaque fois le produit compact Q*bⱼ
        qtprod!(A, view(B, 1:m, j)) #produit vecteur colonne par vecteur colonne de B
    end
    B
end