import Base: getindex, size, setindex!, length, findnz, nnz
export SparseMatrixLIL

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
function SparseMatrixLIL(A :: AbstractMatrix{Tv}) where {Tv}
    SparseMatrixLIL{Tv,Int}(map(v -> sparsevec(A[:,v]), 1:size(A,2)), size(A)...)
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
    for (ji,j) in enumerate(J)
        A.data[j][I] = view(x, :, ji)
    end
    x
end

function nnz(A::SparseMatrixLIL{Tv,Ti}) where {Tv,Ti<:Integer}
    n = 0
    for j = 1:A.n
        n += nnz(A.data[j])
    end
    n
end

function findnz(S::SparseMatrixLIL{Tv,Ti}) where {Tv,Ti}
    numnz = nnz(S)
    I = Vector{Ti}(numnz)
    J = Vector{Ti}(numnz)
    V = Vector{Tv}(numnz)

    count = 1
    for col = 1:S.n
        nzval = S.data[col].nzval
        nzind = S.data[col].nzind
        for u = 1:length(nzind)            
            if nzval[u] != 0
                I[count] = nzind[u]
                J[count] = col
                V[count] = nzval[u]
                count += 1
            end
        end
    end

    count -= 1
    if numnz != count
        deleteat!(I, (count+1):numnz)
        deleteat!(J, (count+1):numnz)
        deleteat!(V, (count+1):numnz)
    end

    return (I, J, V)
end
