export wavescorevote, wavescorematrix

function wavescorevote(m::NullMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    0.0
end
wavescorematrix(m::NullMeasure) = SparseMatrixLIL(spzeros(Float64,size(m.G1,1),size(m.G2,1)))

function wavescorevote(m::ConvexCombMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    m.alpha * wavescorevote(m.S,i,j,adj1,adj2,adjcount1,adjcount2) +
    (1-m.alpha) * wavescorevote(m.T,i,j,adj1,adj2,adjcount1,adjcount2)
end
wavescorematrix(m::ConvexCombMeasure) =
    m.alpha * wavescorematrix(m.S) + (1-m.alpha) * wavescorematrix(m.T)

function wavescorevote(meas::NodeSimMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    0.0
end
wavescorematrix(m::NodeSimMeasure) = copy(m.S)

function wavescorevote(meas::WECMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    meas.S[i,j] .+ meas.S[adj1,adj2]
end
wavescorematrix(m::WECMeasure) = copy(m.S)

function local_cet_ncet(meas::DWECMeasure,
                                   i::Int,j::Int,
                                   adj1::AbstractVector{Int},adj2::AbstractVector{Int})
    m = length(adj1)
    n = length(adj2)
    TC = zeros(Float64,m,n)
    TN = zeros(Float64,m,n)
    @inbounds for u = 1:m, v = 1:n
        Tc,Tn = cet_ncet(meas.G1[i,adj1[u]],meas.G2[j,adj2[v]])
        TC[u,v] = Tc
        TN[u,v] = Tn
    end
    (TC,TN)
end
function wavescorevote(meas::DWECMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    TC,TN = local_cet_ncet(meas, i,j, adj1,adj2)
    activitynorm = let TS = TC+TN; isempty(TS) ? -Inf : maximum(TS) end
    meas.S[i,j] .+ TC .* meas.S[adj1,adj2] ./ activitynorm
end
wavescorematrix(m::DWECMeasure) = copy(m.S)

# Initial matrix = alpha min(G1.expci[i],G2.expci[j])/max(G1.maxdeg,G2.maxdeg) + (1-alpha) S[i,j]
# c1[adj1] += g1.dep[i]; c2[adj2] += g2.dep[j];
# Replace 1st w/ alpha (min(G1.expci[ip]-c1[ip],G2.expci[jp]-c2[jp])/max(G1.maxdeg,G2.maxdeg) + 1)
function wavescorevote(meas::NetalMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    m = length(adj1)
    n = length(adj2)
    R = Matrix{Float64}(m,n)
    for rj = 1:n, ri = 1:m
        u = adj1[ri]
        v = adj2[rj]
        c1 = adjcount1[u] * meas.dep1[u]
        c2 = adjcount2[v] * meas.dep2[v]
        cp1 = c1 - meas.dep1[u]
        cp2 = c2 - meas.dep2[v]
        R[ri,rj] = min(meas.expci1[u]-c1,meas.expci2[v]-c2)/max(meas.maxdeg1,meas.maxdeg2) -
                 min(meas.expci1[u]-cp1,meas.expci2[v]-cp2)/max(meas.maxdeg1,meas.maxdeg2) + 1
    end
    R
end
wavescorematrix(m::NetalMeasure) = expci(m.G1,m.G2)

function wavescorevote(meas::GhostMeasure,
                       i::Int,j::Int,
                       adj1::AbstractVector{Int},adj2::AbstractVector{Int},
                       adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
    m = length(adj1)
    n = length(adj2)
    R = Matrix{Float64}(m,n)
    for rj = 1:n, ri = 1:m
        u = adj1[ri]
        v = adj2[rj]
        R[ri,rj] = score(meas, u, v)
    end
    (score(meas, i, j) .+ R) ./ 2
end
wavescorematrix(m::GhostMeasure) = SparseMatrixLIL(spzeros(Float64,size(m.G1,1),size(m.G2,1)))
