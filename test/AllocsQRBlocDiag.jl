nb_Matrix = 50
m_max = 50
n_max = 30

A_vect, m, n = CreateDiagBlock(nb_Matrix, m_max, n_max)
A = BlockDiagonal(A_vect)
A_aux = deepcopy(A)
QRblocdiag!(A)
QRblocdiag!_a = @allocated begin
    @views QRblocdiag!(A_aux)
end; @show(QRblocdiag!_a)
display(@test QRblocdiag!_a == 0)