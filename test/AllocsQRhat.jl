m  = 100
n = 90
A = rand(m,n)

qrhat_qrH_a = @allocated begin 
    qrhat!(A)
end; #@show(qrhat_qrH_a)
@test qrhat_qrH_a == 0


nb_Matrix = 50
m_max = 50
n_max = 30

A_vect, m, n = CreateDiagBlock(nb_Matrix, m_max, n_max)
Ab = BlockDiagonal(A_vect)
Ab_aux = deepcopy(Ab)
qrhat!(Ab_aux)

qrhat_QRBlocDiag_a = @allocated begin 
    qrhat!(Ab)
end; #@show(qrhat_QRBlocDiag_a)
@test qrhat_QRBlocDiag_a == 0

m = 100
n1 = 50
n2 = 40
n = n1 + n2

Ah = BlockArray(rand(m,n), [m], [n1,n2])

qrhat_QRhcat_a = @allocated begin 
    qrhat!(Ah)
end; #@show(qrhat_QRhcat_a)
#display(@test qrhat_QRhcat_a == 0)