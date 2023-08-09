m = 100
n1 = 50
n2 = 40
n = n1 + n2

A = BlockArray(rand(m,n), [m], [n1,n2])

println("Benchmark for QRhcat! :")
display(@benchmark QRhcat!($A))
println("---------------")