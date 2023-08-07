export QRblocdiag!

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