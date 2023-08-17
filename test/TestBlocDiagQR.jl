
println("-------------------------------")
println("QR factorization for Block-Diagonal matrices")
println("-------------------------------")
    @testset "With the method using Abstract Vectors" begin
        nb_Matrix = 50
        m_max = 50
        n_max = 30

        A_vect, m, n = CreateDiagBlock(nb_Matrix, m_max, n_max)
        R_vect = deepcopy(A_vect)
        QRblocdiag!(R_vect)

        for k in eachindex(A_vect)
            F = qr(A_vect[k])

            m, n = size(R_vect[k])
            Q_H = QRebuild!(R_vect[k])

            @test norm(triu(R_vect[k][1:n,1:n]) - F.R) ≤ 1e-13
            @test norm(Q_H - F.Q) ≤ 1e-13
            @test norm(Q_H*triu(R_vect[k]) - A_vect[k]) ≤ 1e-13
        end
    end

    # Tests made through the use of BlockDiagonals.jl package

    @testset "With the method using Abstract Vectors via BlockDiagonals.jl package" begin
        nb_Matrix = 50
        m_max = 50
        n_max = 30
        v = CreateDiagBlock(nb_Matrix, m_max, n_max)[1]
        A = BlockDiagonal(v)
        A_aux = deepcopy(A)
        v_aux = BlockDiagonals.blocks(A_aux)
        QRblocdiag!(v_aux)
        
        #block by block tests
        for k in eachindex(v)
            F = qr(v[k])

            m, n = size(v_aux[k])
            Q_H = QRebuild!(v_aux[k])

            @test norm(triu(v_aux[k][1:n, 1:n]) - F.R) ≤ 1e-13
            @test norm(Q_H - F.Q) ≤ 1e-13
            @test norm(Q_H*triu(v_aux[k]) - v[k]) ≤ 1e-13
        end
    end

    # Tests using qrhat! with a BlockDiagonal Matrix

    @testset "Testing through the hat function - BlockDiagonal method - Block by block" begin
        nb_Matrix = 50
        m_max = 50
        n_max = 30
        v = CreateDiagBlock(nb_Matrix, m_max, n_max)[1]
        A = BlockDiagonal(v)
        A_aux = deepcopy(A)
        qrhat!(A_aux)
        v_aux = BlockDiagonals.blocks(A_aux)
        
        #block by block tests
        for k in eachindex(v)
            F = qr(v[k])

            m, n = size(v_aux[k])
            Q_H = QRebuild!(v_aux[k])

            @test norm(triu(v_aux[k][1:n, 1:n]) - F.R) ≤ 1e-13
            @test norm(Q_H - F.Q) ≤ 1e-13
            @test norm(Q_H*triu(v_aux[k]) - v[k]) ≤ 1e-13
        end
    end

    @testset "Testing through the hat function - BlockDiagonal method - Full matrix" begin
        bm = BlockDiagonal([rand(3, 3), rand(3, 2)])
        m, n = size(bm)
        am = copy(bm)
        b = rand(m)
        c = copy(b)

        qrhat!(am)
        F = qr(bm)
        R = get_r(am)

        Q = QRebuildBDM!(am)
        qtprod!(am, b)
        @test norm(b - Q'c) ≤ 1e-13
        
        @test norm(F.Q - Q) ≤ 1e-13
        @test norm(F.R - R[1:n, 1:n]) ≤ 1e-13
        @test norm((Q * R) - bm) ≤ 1e-13
    end