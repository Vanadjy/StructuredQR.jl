m  = 100
n = 90
A = rand(m,n)

println("Benchmark for qrhat! - qrH! :")
display(@benchmark qrhat!($A))
println("---------------")

nb_Matrix = 50
m_max = 50
n_max = 30

A_vect, m, n = CreateDiagBlock(nb_Matrix, m_max, n_max)
Ab = BlockDiagonal(A_vect)

println("Benchmark for qrhat! - QRblocdiag! :")
display(@benchmark qrhat!($Ab))
println("---------------")

m = 100
n1 = 50
n2 = 40
n = n1 + n2

Ah = BlockArray(rand(m,n), [m], [n1,n2])

println("Benchmark for qrhat! - QRhcat! :")
display(@benchmark qrhat!($Ah))
println("---------------")