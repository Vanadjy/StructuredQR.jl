export qrhat!

qrhat!(A::AbstractMatrix) = qrH!(A) #Dense matrix

qrhat!(A::AbstractSparseMatrix) = LinearAlgebra.qr!(A) #Sparse matrix

qrhat!(A::BlockDiagonal{T, Matrix{T}}) where T = QRblocdiag!(blocks(A)) #BlockDiagonal matrix

qrhat!(A::HcatMatrix{Matrix{T}, Matrix{T}}) where T = QRhcat!(A.leftm, A.rightm) #Horizontally concatenated matrix