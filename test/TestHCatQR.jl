## Tests of QR with horizontal concatenation

println("-------------------------------")
println("QR factorization for Horizontally concatenated matrices")
println("-------------------------------")
## Overdetermined matrices ##
@testset "Overdetermined matrices - Method with two matrices" begin
    m = 6
    n1 = 2
    n2 = 3
    n = n1 + n2

    A = rand(m, n)
    A1, A2 = splitH(A, n1)

    B = deepcopy(A)
    B1 = similar(B)
    B2 = similar(B)
    F = qr(A)

    R_H = deepcopy(A)
    qrH!(R_H)
    Q_H = QRebuild!(R_H)
    qmul!(B1, R_H, B)
    @test norm(Q_H * B - B1) <= 1e-13 #tests if the multiplication QB is correct
    qtmul!(B2, R_H, B)
    @test norm(Q_H'B - B2) <= 1e-13 #tests if the multiplication Q*B is correct

    QRhcat!(A1, A2)
    R = RRebuildHcat(A1, A2)
    Q = QRebuildHcat(A1, A2)
    @test norm(F.R - R[1:n, 1:n]) <= 1e-14 #test of unicity of QR decomposition
    @test norm(F.Q - Q) <= 1e-13
    @test norm(A - Q * R) <= 1e-13 #tests if the decomposition is correct
end

## Square matrices ##
@testset "Square matrices - Method with two matrices" begin
    m = 6
    n1 = 3
    n2 = 3
    n = n1 + n2

    A = rand(m, n)
    A1, A2 = splitH(A, n1)
    B = deepcopy(A)
    B1 = similar(B)
    B2 = similar(B)
    F = qr(A)

    R_H = deepcopy(A)
    qrhat!(R_H)
    Q_H = QRebuild!(R_H)
    qmul!(B1, R_H, B)
    @test norm(Q_H * B - B1) <= 1e-13 #tests if the multiplication QB is correct
    qtmul!(B2, R_H, B)
    @test norm(Q_H'B - B2) <= 1e-13 #tests if the multiplication Q*B is correct

    QRhcat!(A1, A2)
    R = RRebuildHcat(A1, A2)
    Q = QRebuildHcat(A1, A2)
    @test norm(F.R - R[1:n, 1:n]) <= 1e-13 #test of unicity of QR decomposition
    @test norm(F.Q - Q) <= 1e-13
    @test norm(A - Q * R) <= 1e-13 #tests if the decomposition is correct
end

# Tests using qrhat!
@testset "Testing through the hat function - BlockArray method" begin
    m = 6
    n1 = 2
    n2 = 3
    n = n1 + n2

    A1, A2 = rand(m, n1), rand(m, n2)
    R_H = BlockArray(rand(m, n), [m], [n1, n2])
    A = Array(R_H)
    A_aux = deepcopy(A)
    B = deepcopy(A)
    B1 = similar(B)
    B2 = similar(B)
    F = qr(A)

    qrhat!(R_H)
    Q_H = QRebuild!(R_H)
    qmul!(B1, R_H, B)
    @test norm(Q_H * B - B1) <= 1e-13 #tests if the multiplication QB is correct
    qtmul!(B2, R_H, B)
    @test norm(Q_H'B - B2) <= 1e-13 #tests if the multiplication Q*B is correct

    R = RRebuildHcat(R_H[Block(1, 1)], R_H[Block(1, 2)])
    Q = QRebuildHcat(R_H[Block(1, 1)], R_H[Block(1, 2)])
    @test norm(F.R - R[1:n, 1:n]) <= 1e-14 #test of unicity of QR decomposition
    @test norm(F.Q - Q) <= 1e-13
    @test norm(A - Q * R) <= 1e-13 #tests if the decomposition is correct
end