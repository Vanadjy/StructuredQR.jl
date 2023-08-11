m = 100
n1 = 50
n2 = 40
n = n1 + n2

A = rand(m, n)
b = rand(m)
B = rand(m,n)
R = copy(A)
qrH!(R)

qtprod!_a = @allocated begin
    qtprod!(R,b)
end; #@show(qtprod!_a)
@test qtprod!_a == 0

qprod!_a = @allocated begin
    qprod!(R,b)
end; #@show(qprod!_a)
@test qprod!_a == 0

qtmul!_a = @allocated begin
    qtmul!(R,B)
end; #@show(qtmul!_a)
@test qtmul!_a == 0

qmul!_a = @allocated begin
    qtmul!(R,B)
end; #@show(qmul!_a)
@test qmul!_a == 0