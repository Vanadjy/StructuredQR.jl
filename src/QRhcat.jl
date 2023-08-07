export QRhcat!

function QRhcat!(A1::AbstractMatrix, A2::AbstractMatrix)
    """
    QRhcat!(A1::AbstractMatrix, A2::AbstractMatrix)

    Computes the QR-factorisation of an overdetermined full-rank (m × n) matrix A presented as a horizontally conctenated one, i.e. A = [A₁ A₂] where A₁ and A₂ are overdetermined full-rank matrices of dimension (m × n₁) and (m × n₂) respectively such that : n₁ + n₂ = n.

    Computes : A = QR

    Where :
        - Q is an unitary matrix (it means Q*Q = I)
        - R is upper triangular
    
    As A is seen as horizontally concatenated, we will first compute A₁ = Q₁R₁ to then calculate Q₁ᵀA₂ and compute afterward the dense QR factorization on the below block of Q₁ᵀA₂ (of dimension (m-n₁ × n₂)

    #### Input arguments

    * `A1` : a full-rank overdetermined matrix of dimension (m × n₁)
    * `A2` : a full-rank overdetermined matrix of dimension (m × n₂)

    #### Output arguments

    * `A1` : a full-rank overdetermined matrix of dimension (m × n₁), containing the coefficients of Q₁ and R₁
    * `A2` : a full-rank overdetermined matrix of dimension (m × n₂), modified such that A = QR
    """
    n1 = size(A1,2) #nombre de colonne du premier bloc
    m2, n2 = size(A2)

    qrH!(A1) #compute QR Householder on the first block of A : A1

    qtmul!(A1, A2) #compute Q' *A2

    Q1_orth_A2 = view(A2,n1+1:m2,1:n2) #Q1_orth_A2 = A[n1+1:end, n1+1:end]
    qrH!(Q1_orth_A2) #compute QR Householder on Q1_orth_A2 which is the second block of Q1*A2

    return A1, A2
end