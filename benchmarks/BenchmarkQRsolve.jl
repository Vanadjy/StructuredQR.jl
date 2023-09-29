n1 = 5
n2 = 3
m = 10
n = n1 + n2
A = rand(m, n)
b = rand(m)
res = deepcopy(b)

R = get_r(A)

println("Benchmark for StructuredQR.rdiv! :")
display(@benchmark StructuredQR.rdiv!($res, $R, $b))
println("---------------")

println("Benchmark for qrsolve! :")
display(@benchmark qrsolve!($res, $R, $b))
println("---------------")