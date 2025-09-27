m = 100
n1 = 50
n2 = 40
n = n1 + n2

A = BlockArray(rand(m,n), [m], [n1,n2])

QRhcat!_a = @allocated begin
    QRhcat!(A)
end; #@show(QRhcat!_a)
@test QRhcat!_a == 0