export my_sign, HouseholderReflection, splitH, CreateDiagBlock

function HouseholderReflection(u::AbstractVector)
    """
    HouseholderReflection(u::AbstractVector)

    Builds the Householder reflector associated with the vector u

    #### Input arguments

    * `u` : a vector of size n

    #### Output arguments

    * `H` where H = I - 2(uu*/u*u) the (n × n) Householder reflector, a rank 1 matrix
    """
    δ = u'u
    return (I - 2*u*u'/δ)
end

function splitH(A::AbstractMatrix, n1::Int)
    return A[1:end, 1:n1], A[1:end, n1+1:end]
end

function CreateDiagBlock(nb_Matrix::Int, m_max::Int, n_max::Int)
    A_vect = Matrix{Float64}[]
    m_vect = rand(1:m_max, nb_Matrix)
    n_vect = rand(1:n_max, nb_Matrix)
    for k = 1:nb_Matrix
        if m_vect[k]<n_vect[k] #ensures to always being in overdetermined cases
            (m_vect[k],n_vect[k]) = (n_vect[k],m_vect[k])
        end
        push!(A_vect, rand(m_vect[k], n_vect[k]))
    end
    m = sum(m_vect)
    n = sum(n_vect)
    return A_vect, m , n 
end