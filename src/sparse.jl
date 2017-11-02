import Base: getindex, size, setindex!, length, findnz

# Store sparse matrix as as vector of sparse vectors (list of lists)
# This is to enable more efficient updates (i.e. setindex!)
# Each vector represents a COLUMN on the matrix
immutable SparseMatrixLIL{Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}
    data :: Vector{SparseVector{Tv,Ti}}
    m :: Int
    n :: Int
end
function SparseMatrixLIL(A :: SparseMatrixCSC{Tv,Ti}) where {Tv,Ti<:Integer}
    SparseMatrixLIL{Tv,Ti}(map(v -> A[:,v], 1:size(A,2)), size(A)...)
end
function SparseMatrixLIL(A :: AbstractMatrix)
    SparseMatrixLIL(map(v -> sparsevec(A[:,v]), 1:size(A,2)), size(A)...)
end
size(A::SparseMatrixLIL) = (A.m,A.n)
size(A::SparseMatrixLIL,d::Integer) = (A.m,A.n)[d]
length(A::SparseMatrixLIL) = A.m * A.n
getindex(A::SparseMatrixLIL, i::Integer, j::Integer) = A.data[j][i]
getindex(A::SparseMatrixLIL, I::AbstractVector, j::Integer) = SparseMatrixLIL([A.data[j][I]])
function getindex(A::SparseMatrixLIL{Tv,Ti}, I::AbstractVector, J::AbstractVector) where {Tv,Ti<:Integer}
    data = Vector{SparseVector{Tv,Ti}}(length(J))
    for (ji,j) in enumerate(J)
        data[ji] = A.data[j][I]
    end
    SparseMatrixLIL{Tv,Ti}(data, length(I), length(J))
end
function setindex!(A::SparseMatrixLIL, x, i::Integer, j::Integer)
    A.data[j][i] = x
end
function setindex!(A::SparseMatrixLIL, x::AbstractVector, I::AbstractVector, j::Integer)
    A.data[j][I] = x
end
function setindex!(A::SparseMatrixLIL{Tv,Ti}, x::AbstractMatrix, I::AbstractVector, J::AbstractVector) where {Tv,Ti<:Integer}
    data = Vector{SparseVector{Tv,Ti}}(length(J))
    for (ji,j) in enumerate(J)
        A.data[j][I] = view(x, :, ji)
    end
    x
end
function findnz(A::SparseMatrixLIL{Tv,Ti}) where {Tv,Ti<:Integer}
    I = Ti[]
    J = Ti[]
    V = Tv[]
    for j = 1:A.n
        lst = A.data[j].nzval
        for (u,i) in enumerate(A.data[j].nzind)
            push!(I,i)
            push!(J,j)
            push!(V,lst[u])
        end
    end
    (I,J,V)
end
