@testset "rdiv! - dense matrices" begin
    n1 = 5
    n2 = 3
    m = 10
    n = n1 + n2
    A = rand(m, n)
    b = rand(m)
    c = copy(b)
    B = copy(A)

    qrhat!(A)
    F = qr(B)
    R = get_r(A)

    StructuredQR.rdiv!(R, b)
    y = R\c
    @test norm(y - b[1:n]) ≤ 1e-12
end

@testset "rdiv! - Horizontally concatenated matrices" begin
    n1 = 5
    n2 = 3
    m = 10
    n = n1 + n2
    A = BlockArray(rand(m,n), [m], [n1,n2])
    b = rand(m)
    c = copy(b)
    B = copy(A)

    qrhat!(A)
    F = qr(B)
    R = get_r(A)

    StructuredQR.rdiv!(R, b)
    y = R\c
    @test norm(y - b[1:n]) ≤ 1e-12
end

@testset "rdiv! - BlockDiagonal matrices" begin
    nb_Matrix = 50
    m_max = 50
    n_max = 30
    
    A_vect = CreateDiagBlock(nb_Matrix, m_max, n_max)[1]
    A = BlockDiagonal(A_vect)
    B = copy(A)

    qrhat!(A)
    vA = BlockDiagonals.blocks(A)
    vB = BlockDiagonals.blocks(B)

    for k in eachindex(vA)
        m, n = size(vA[k])
        F = qr(vB[k])
        R = get_r(vA[k])
        b = rand(m)
        c = copy(b)
    
        StructuredQR.rdiv!(R, b)
        y = R\c
        @test norm(y - b[1:n]) ≤ 1e-12
    end
end

@testset "qrsolve! - BlockDiagonal matrices - Full matrix" begin
    bm = BlockDiagonal([rand(3, 3), rand(3, 2)])
    m, n = size(bm)
    b = rand(m)
    c = copy(b)
    am = copy(bm)

    qrhat!(am)
    F = qr(bm)
    R = get_r(am)

    StructuredQR.rdiv!(am, b)
    y = R\c
    @test norm(b[1:n] - y) ≤ 1e-12
end

@testset "qrsolve! - dense matrices" begin
    n1 = 5
    n2 = 3
    m = 10
    n = n1 + n2
    A = rand(m, n)
    b = rand(m)
    c = copy(b)
    B = copy(A)
    
    F = qr(B)
    
    qrsolve!(A, b)
    y = F\c
    @test norm(y - b[1:n]) ≤ 1e-12
end

@testset "qrsolve! - horizontally concatenated matrices" begin
    n1 = 5
    n2 = 3
    m = 10
    n = n1 + n2
    A = BlockArray(rand(m,n), [m], [n1,n2])
    b = rand(m)
    c = copy(b)
    B = copy(A)
    
    F = qr(B)
    
    qrsolve!(A, b)
    y = F\c
    @test norm(y - b[1:n]) ≤ 1e-12
end

@testset "qrsolve! - BlockDiagonal matrices - block by block" begin
    nb_Matrix = 50
    m_max = 50
    n_max = 30
    
    A_vect = CreateDiagBlock(nb_Matrix, m_max, n_max)[1]
    A = BlockDiagonal(A_vect)
    B = copy(A)

    vA = BlockDiagonals.blocks(A)
    vB = BlockDiagonals.blocks(B)

    for k in eachindex(vA)
        m, n = size(vA[k])
        F = qr(vB[k])

        b = rand(m)
        c = copy(b)
    
        qrsolve!(vA[k], b)
        y = F\c
        @test norm(y - b[1:n]) ≤ 1e-12
    end
end

@testset "qrsolve! - BlockDiagonal matrices - Full matrix" begin
    bm = BlockDiagonal([rand(3, 3), rand(3, 2)])
    m, n = size(bm)
    b = rand(m)
    c = copy(b)
    am = copy(bm)
    qrsolve!(am, b)
    F = qr(bm)
    y = F\c
    @test norm(b[1:n] - y) ≤ 1e-12
end