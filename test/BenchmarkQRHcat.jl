m = 100
n1 = 50
n2 = 40
n = n1 + n2

A = rand(m, n)
b = rand(m)
B = rand(m,n)
A1, A2 = splitH(A, n1)

println("Benchmark for QRhcat! :")
display(@benchmark QRhcat!($A1, $A2))
println("---------------")