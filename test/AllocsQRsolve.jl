n1 = 5
n2 = 3
m = 10
n = n1 + n2
A = rand(m, n)
b = rand(m)
res = deepcopy(b)

R = get_r(A)

a_rdiv! = @allocated begin
    StructuredQR.rdiv!(res, R, b)
end; #show(a_rdiv!)
@test a_rdiv! == 0

a_qrsolve! = @allocated begin
    qrsolve!(res, A, b)
end; #show(a_qrsolve!)
@test a_qrsolve! == 0