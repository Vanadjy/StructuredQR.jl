export QRblocdiag!, get_block

using BlockDiagonals

function get_block(A::BlockDiagonal{T, Matrix{T}}, k::Int) where T
    #Idea : return the kᵗʰ block of A without allocating a new vector through BlockDiagonals.blocks
    index_row = 0
    index_col = 0
    for l=1:(k-1)
        index_row += BlockDiagonals.blocksize(A,l)[1]
        index_col += BlockDiagonals.blocksize(A,l)[2]
    end
    (mₖ, nₖ) = BlockDiagonals.blocksize(A,k)
    return view(A,index_row+1:index_row+mₖ, index_col+1:index_col+nₖ)
end

function QRblocdiag!(A_vect::AbstractVector)
    """
    QRblocdiag!(A_vect::AbstractVector)

    Computes the QR factorization of an overdetermined Block-Diagonal matrix A through all of its full-rank overdetermined diagonal blocks A₁, A₂, ... Aᵣ as A = QR where Q = diag(Q₁, Q₂, ..., Qᵣ) and R = diag(R₁, R₂, ..., Rᵣ) such that :
    
    A₁ = Q₁R₁
    A₂ = Q₂R₂
    ...
    Aᵣ = QᵣRᵣ

    Each Aᵢ block stores the QR information just like qrH! does : the upper triangular stores the coefficients of Rᵢ and the other coefficients of Aᵢ store the vectors required to rebuild Qᵢ

    Computes : A = QR

    Where : 
        - Q is an unitary Block-Diagonal matrix
        - R is an upper triangular Block-Diagonal matrix
    
    #### Input arguments

    * `A_vect` : an AbstractVector of size r containing all the blocks of a Block-Diagonal matrix;

    #### Output arguments

    * `A_vect` : the same vector where all the matrices have been QR-factorized;
    """

    for A in A_vect
        qrH!(A)
    end
end

function QRblocdiag!(A::BlockDiagonal{T, Matrix{T}}) where T
    """
    QRblocdiag!(A::BlockDiagonal{T, Matrix{T}})

    Computes the QR factorization of an overdetermined Block-Diagonal matrix A through all of its full-rank overdetermined diagonal blocks A₁, A₂, ... Aᵣ as A = QR where Q = diag(Q₁, Q₂, ..., Qᵣ) and R = diag(R₁, R₂, ..., Rᵣ) such that :
    
    A₁ = Q₁R₁
    A₂ = Q₂R₂
    ...
    Aᵣ = QᵣRᵣ

    Each Aᵢ block stores the QR information just like qrH! does : the upper triangular stores the coefficients of Rᵢ and the other coefficients of Aᵢ store the vectors required to rebuild Qᵢ

    Computes : A = QR

    Where : 
        - Q is an unitary Block-Diagonal matrix
        - R is an upper triangular Block-Diagonal matrix
    
    #### Input arguments

    * `A` : a BlockDiagonal{T, Matrix{T}} defined by its r blocks;

    #### Output arguments

    * `A` : the same BlockDiagonal Matrix whose r blocks now contain its own QR factorization information;
    """
    for j in 1:(nblocks(A))
        qrH!(BlockDiagonals.blocks(A)[j])
    end
    A
end