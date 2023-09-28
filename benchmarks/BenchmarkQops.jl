m = 100
n1 = 50
n2 = 40
n = n1 + n2

A = rand(m, n)
b = rand(m)
res = similar(b)
B = rand(m,n)
RES = similar(B)

R = copy(A)
qrH!(R)

## qtprod! ##

println("Benchmark for qtprod! :")
display(@benchmark qtprod!($R, $b))
println("---------------")

println("Benchmark for qtprod! :")
display(@benchmark qtprod!($res, $R, $b))
println("---------------")

## qprod! ##

println("Benchmark for qprod! :")
display(@benchmark qprod!($R, $b))
println("---------------")

println("Benchmark for qprod! :")
display(@benchmark qprod!($res, $R, $b))
println("---------------")

## qtmul! ##

println("Benchmark for qtmul! :")
display(@benchmark qtmul!($R, $B))
println("---------------")

println("Benchmark for qtmul! :")
display(@benchmark qtmul!($RES, $R, $B))
println("---------------")

## qmul! ##

println("Benchmark for qmul! :")
display(@benchmark qmul!($R, $B))
println("---------------")

println("Benchmark for qmul! :")
display(@benchmark qmul!($RES, $R, $B))
println("---------------")