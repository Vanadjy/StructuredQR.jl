export QRebuild!, QRebuildHcat, RRebuildHcat

function QRebuild!(A::AbstractMatrix; Q=I)
    """
    QRebuild!(A::AbstractMatrix; Q=I)

    Rebuilds the unitary matrix Q from the information stored in the matrix A that has been transformed by the function Householder_Compact! so that A = QR
        
    Computes : Q so that Q is unitary (i.e. Q*Q = I)

    As the vectors stored in A are scaled, we can just get build the vector [1 ; A[j+1:end,j]] since, as it has been explained in the doc string of the function Householder_Compact!, the first element is supposed to be 1 and the others are scaled with respect to this. 

    We build the matrix Hⱼ which is the Reflexion of Householder of dimension (m-j+1)x(m-j+1) and Q is computed one step after the other according to the following rule :

    Q ← Q × Qⱼ where Qⱼ = [Iⱼ 0]
                          [0 Hⱼ]

    #### Input arguments

    * `A`: a full rank matrix of dimension m × n;

    #### Keyword arguments

    * `Q=I`: the identity matrix of dimension m × m , from which the rebuilding begins;
                      
    #### Output arguments

    * `Q`: the unitary matrix that fits with the QR decomposition of A

    NB : the purpose of this function is essentially to carry out tests to check if our QR decomposition verifies its uniqueness property (depending on the sign of the diagonal terms of R)
    """
    m, n = size(A)
    j = 1
    while (j <= n) & (j < m)
        u_j = [1 ; A[j+1:end,j]]
        Hj = HouseholderReflection(u_j)
        Qj = [I zeros(j-1,m-j+1) ;
            zeros(m-j+1,j-1) Hj]
        Q = Q*Qj
        j += 1
    end
    Q
end

function QRebuildHcat(A1::AbstractMatrix, A2::AbstractMatrix) #A1 et A2 sont les blocs modifiés dans QR_concat_h_compact!
    """
    QRebuildHcat(A1::AbstractMatrix, A2::AbstractMatrix)

    Rebuilds the unitary matrix Q from the information stored in the matrices A₁ and A₂ that have been transformed by the function qrH! so that A = [A₁ A₂] ; A = QR
        
    Computes : Q so that Q is unitary (i.e. Q*Q = I)

    As the vectors stored in A are scaled, we can just get build the vector [1 ; A[j+1:end,j]] since, as it has been explained in the doc string of the function Householder_Compact!, the first element is supposed to be 1 and the others are scaled with respect to this. 

    We build the matrix Hⱼ which is the Reflexion of Householder of dimension (m-j+1)x(m-j+1) and Q is computed one step after the other according to the following rule :

    Q ← Q × Qⱼ where Qⱼ = [Iⱼ 0]
                          [0 Hⱼ]

    #### Input arguments

    * `A`: a full rank matrix of dimension m × n;

    #### Keyword arguments

    * `Q=I`: the identity matrix of dimension m × m , from which the rebuilding begins;
                      
    #### Output arguments

    * `Q`: the unitary matrix that fits with the QR decomposition of A

    NB : the purpose of this function is essentially to carry out tests to check if our QR decomposition verifies its uniqueness property (depending on the sign of the diagonal terms of R)
    """
    
    m, n1 = size(A1)
    Q1 = QRebuild!(A1)
    Q_prime = QRebuild!(A2[n1+1:end, 1:end])
    return Q1*[I zeros(n1, m-n1); zeros(m-n1, n1) Q_prime]
end

function RRebuildHcat(A1::AbstractMatrix, A2::AbstractMatrix)
    n1 = size(A1, 2)
    R1 = triu(A1)
    R_prime = triu(A2[n1+1:end, 1:end])
    R2 = [A2[1:n1, 1:end] ; R_prime]
    return [R1 R2]
end

function QRebuildBDM!(A::BlockDiagonal{T, Matrix{T}}; Q_vect = Matrix{Float64}[]) where T
    v = BlockDiagonals.blocks(A)
    for M in v
        push!(Q_vect, QRebuild!(M))
    end
    return BlockDiagonal(Q_vect)
end