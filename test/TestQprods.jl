# Tests of Q multiplications where Q is from a QR-factorization

println("-------------------------------")
println("Q multiplications from dense matrices")
println("-------------------------------")
## Overdetermined matrices ##
@testset "Overdetermined Matrices - Method (res, A, b)" begin
    m, n = 10, 8
    A = rand(m,n)
    b = rand(m)
    b1 = similar(b)
    b2 = similar(b)

    B = copy(A)
    B1 = similar(B)
    B2 = similar(B1)

    R_H = deepcopy(A)
    qrhat!(R_H)
    Q_H = QRebuild!(R_H)

    qtprod!(b1, R_H, b)
    @test norm((Q_H')*b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(b2, R_H, b)
    @test norm((Q_H)*b - b2) <= 1e-13
    qmul!(B1, R_H, B)
    @test norm(Q_H*B - B1) <= 1e-13
    qtmul!(B2, R_H, B)
    @test norm(Q_H'B - B2) <= 1e-13
end

@testset "Overdetermined Matrices - Method (A, b)" begin
    m, n = 10, 8
    A = rand(m,n)
    b = rand(m)
    b1 = copy(b)
    b2 = copy(b)

    B = copy(A)
    B1 = copy(B)
    B2 = copy(B1)

    R_H = deepcopy(A)
    qrhat!(R_H)
    Q_H = QRebuild!(R_H)

    qtprod!(R_H, b1)
    @test norm((Q_H')*b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(R_H, b2)
    @test norm((Q_H)*b - b2) <= 1e-13
    qmul!(R_H, B1)
    @test norm(Q_H*B - B1) <= 1e-13
    qtmul!(R_H, B2)
    @test norm(Q_H'B - B2) <= 1e-13
end

## Square matrices ##
@testset "Square Matrices - Method (res, A, b)" begin
    m, n = 10, 10
    A = rand(m,n)
    b = rand(m)
    b1 = deepcopy(b)
    b2 = deepcopy(b)

    B = copy(A)
    B1 = copy(B)
    B2 = copy(B1)

    R_H = deepcopy(A)
    qrhat!(R_H)
    Q_H = QRebuild!(R_H)

    qtprod!(b1, R_H, b)
    @test norm((Q_H')*b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(b2, R_H, b)
    @test norm((Q_H)*b - b2) <= 1e-13
    qmul!(B1, R_H, B)
    @test norm(Q_H*B - B1) <= 1e-13
    qtmul!(B2, R_H, B)
    @test norm(Q_H'B - B2) <= 1e-13
end

@testset "Square Matrices - Method (A, b)" begin
    m, n = 10, 10
    A = rand(m,n)
    b = rand(m)
    b1 = deepcopy(b)
    b2 = deepcopy(b)

    B = copy(A)
    B1 = copy(B)
    B2 = copy(B1)

    R_H = deepcopy(A)
    qrhat!(R_H)
    Q_H = QRebuild!(R_H)

    qtprod!(R_H, b1)
    @test norm((Q_H')*b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(R_H, b2)
    @test norm((Q_H)*b - b2) <= 1e-13
    qmul!(R_H, B1)
    @test norm(Q_H*B - B1) <= 1e-13
    qtmul!(R_H, B2)
    @test norm(Q_H'B - B2) <= 1e-13
end

println("-------------------------------")
println("Q multiplications from dense matrices - LinearOperator")
println("-------------------------------")

## Overdetermined matrices ##
@testset "Overdetermined Matrices" begin
    m, n = 10, 8
    A = rand(m,n)
    b = rand(m)
    b1 = deepcopy(b)
    b2 = deepcopy(b)
    R_H = deepcopy(A)

    qrhat!(R_H)
    op = q_lop(R_H)

    qtprod!(b1, R_H, b)
    @test norm(op'b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(b2, R_H, b)
    @test norm(op*b - b2) <= 1e-13
end

## Square matrices ##
@testset "Square Matrices" begin
    m, n = 10, 10
    A = rand(m,n)
    b = rand(m)
    b1 = deepcopy(b)
    b2 = deepcopy(b)
    R_H = deepcopy(A)

    qrhat!(R_H)
    op = q_lop(R_H)

    qtprod!(b1, R_H, b)
    @test norm(op'b - b1) <= 1e-13 #tests if the multiplication is correct
    qprod!(b2, R_H, b)
    @test norm(op*b - b2) <= 1e-13
end


println("-------------------------------")
println("Q multiplications from BlockDiagonal matrices")
println("-------------------------------")

@testset "Testing through the hat function - BlockDiagonal method - Full matrix" begin
    bm = BlockDiagonal([rand(3, 3), rand(3, 2)])
    m, n = size(bm)
    am = copy(bm)
    b = rand(m)
    b1 = similar(b)
    b2 = similar(b)

    qrhat!(am)
    Q = QRebuildBDM!(am)
    qtprod!(b1, am, b)
    @test norm(b1 - Q'b) ≤ 1e-13
    qprod!(b2, am, b)
    @test norm(b2 - Q*b) ≤ 1e-13
end

println("-------------------------------")
println("Q multiplications from BlockDiagonal matrices - LinearOperator")
println("-------------------------------")

@testset "Overdetermined Matrices" begin
    bm = BlockDiagonal([rand(3, 3), rand(3, 2)])
    m, n = size(bm)
    am = copy(bm)
    b = rand(m)
    b1 = similar(b)
    b2 = similar(b)

    qrhat!(am)
    op = q_lop(am)
    qtprod!(b1, am, b)
    @test norm(b1 - op'b) ≤ 1e-13
    qprod!(b2, am, b)
    @test norm(b2 - op*b) ≤ 1e-13
end