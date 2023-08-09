export qrhat!

qrhat!(A::AbstractMatrix) = qrH!(A) #Dense matrix

qrhat!(A::AbstractSparseMatrix) = LinearAlgebra.qr!(A) #Sparse matrix

qrhat!(A::BlockDiagonal{T, Matrix{T}}) where T = QRblocdiag!(A) #BlockDiagonal matrix

qrhat!(A::AbstractBlockMatrix{T}) where T = QRhcat!(A) #Horizontally concatenated matrix