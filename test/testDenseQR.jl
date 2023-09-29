## Tests of dense Householder QR factorization

println("-------------------------------")
println("QR factorization for dense matrices")
println("-------------------------------")
    ## Overdetermined matrices ##
    @testset "Overdetermined Matrices" begin
        m, n = 10, 8
        A = rand(m,n)
        R_H = deepcopy(A)
        qrH!(R_H)
        F = qr(A)
        
        Q_H = QRebuild!(R_H)

        @test norm(Q_H'Q_H - I) <= 1e-13 #tests that Q_H si unitary
        @test norm(F.Q - Q_H) <= 1e-13 #tests of unicity of QR decomposition
        @test norm(F.R - triu(R_H[1:n,1:n])) <= 1e-13
        @test norm(Q_H*triu(R_H) - A) <= 1e-13 #tests if the QR decomposition is correct
    end

    ## Square matrices ##
    @testset "Square Matrices" begin
        m, n = 10, 10
        A = rand(m,n)
        R_H = deepcopy(A)
        qrH!(R_H)
        F = qr(A)
        
        Q_H = QRebuild!(R_H)

        @test norm(Q_H'Q_H - I) <= 1e-13 #tests that Q_H si unitary
        @test norm(F.Q - Q_H) <= 1e-13 #tests of unicity of QR decomposition
        @test norm(F.R - triu(R_H[1:n,1:n])) <= 1e-13
        @test norm(Q_H*triu(R_H) - A) <= 1e-13 #tests if the QR decomposition is correct
    end

    # Tests using qrhat!

    @testset "Testing through the hat function" begin
        m, n = 10, 8
        A = rand(m,n)
        R_H = deepcopy(A)
        qrhat!(R_H)
        F = qr(A)
        
        Q_H = QRebuild!(R_H)

        @test norm(Q_H'Q_H - I) <= 1e-13 #tests that Q_H si unitary
        @test norm(F.Q - Q_H) <= 1e-13 #tests of unicity of QR decomposition
        @test norm(F.R - triu(R_H[1:n,1:n])) <= 1e-13
        @test norm(Q_H*triu(R_H) - A) <= 1e-13 #tests if the QR decomposition is correct
    end