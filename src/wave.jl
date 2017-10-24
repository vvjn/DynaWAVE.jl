using DataStructures
using UnicodePlots

export wave, dynawave, waveshuffle, WaveModel, align!

type WaveModel{E,T<:AbstractMatrix,M<:NetalignMeasure}
    G1 :: SparseMatrixCSC{E,Int}
    G2 :: SparseMatrixCSC{E,Int}
    D :: T # vote matrix
    f :: Vector{Int} # mapping
    phi :: M # objective value
    seeds :: Array{Int} # initial seeds k x 2 int array
    objval :: Float64 # objective value
    objval_edge :: Float64
    objval_node :: Float64
    history :: Vector{Float64} # approx. obj. val. at each iteration

    function WaveModel{E,T,M}(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                       D::AbstractMatrix,
                       ecmeasure::NetalignMeasure,
                       seeds) where {E,T,M}
        new(G1,G2, D,
            zeros(Int,size(G1,1)),
            ecmeasure,seeds,0.0,0.0,0.0,Float64[])
    end
end
function WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                   ecmeasure::NetalignMeasure, seeds = zeros(Int,(0,0)))
    D = wavescorematrix(ecmeasure)
    WaveModel{eltype(G1),typeof(D),typeof(ecmeasure)}(G1,G2,D,ecmeasure,seeds)
end

function initializewave!(M::WaveModel;approxobjval=false)
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    # adjcount[i] is the num. of times that a neighbor of i is aligned
    adjcount1 = zeros(Int,n1)
    adjcount2 = zeros(Int,n2)
    M.objval = 0.0
    L1 = Set{Int}()
    L2 = Set{Int}() # nodes already aligned
    for k = 1:size(M.seeds,1)
        i = M.seeds[k,1]
        j = M.seeds[k,2]
        adj1 = adjnodes(M.G1,i)
        adj2 = adjnodes(M.G2,j)
        M.f[i] = j
        push!(L1,i)
        push!(L2,j)
        adjcount1[adj1] += 1
        adjcount2[adj2] += 1
        # add vote to adjacent nodes of both networks
        Ddiff = wavescorevote(M.phi, i,j, adj1,adj2, adjcount1,adjcount2)
        M.D[adj1,adj2] += Ddiff

        if approxobjval
            M.objval += sum(Ddiff[findin(adj1, L1), findin(adj2, L2)])
            push!(M.history, M.objval)
        end
    end
    (L1,L2,adjcount1,adjcount2)
end

function buildpriorityqueue(M::WaveModel,L1::Set{Int},L2::Set{Int})
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    Q = PriorityQueue{Tuple{Int,Int},Float64}()
    if issparse(M.D)
        I,J,V = findnz(M.D)
        for u = 1:length(V)
            i = I[u]
            j = J[u]
            if !(i in L1 || j in L2)
                enqueue!(Q, (i,j), -V[u])
            end
        end
    else
        for i = 1:n1, j = 1:n2
            if !(i in L1 || j in L2)
                enqueue!(Q, (i,j), -M.D[i,j])
            end
        end
    end
    Q
end

function updatepriorityqueue!(M::WaveModel,Q::PriorityQueue,L1::Set{Int,},L2::Set{Int},
                              i::Int,j::Int,
                              adj1::AbstractVector{Int},adj2::AbstractVector{Int})
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    # remove cross node pairs containing i or j
    if issparse(M.D)
        for ip = setdiff(1:n1,L1)
            if haskey(Q,(ip,j)) dequeue!(Q,(ip,j)) end
        end
        for jp = setdiff(1:n2,L2)
            if haskey(Q,(i,jp)) dequeue!(Q,(i,jp)) end
        end
    else
        for ip = setdiff(1:n1,L1); dequeue!(Q,(ip,j)); end
        for jp = setdiff(1:n2,L2); dequeue!(Q,(i,jp)); end
    end
    # update priority queue with the new votes
    if issparse(M.D)
        Ip = setdiff(adj1,L1)
        Jp = setdiff(adj2,L2)
        Vp = -M.D[Ip,Jp]
        for u = 1:length(Ip), v = 1:length(Jp)
            Q[(Ip[u],Jp[v])] = Vp[u,v]
        end
    else
        for ip in setdiff(adj1,L1), jp in setdiff(adj2,L2)
            Q[(ip,jp)] = -M.D[ip,jp]
        end
    end
end

function align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-size(M.seeds,1)),
                approxobjval::Bool=false)
    # intialize alignment using the seeds aligned pairs
    L1,L2,adjcount1,adjcount2 = initializewave!(M,approxobjval=approxobjval)

    # initialize priority queue
    println("Building priority queue")
    Q = buildpriorityqueue(M,L1,L2)

    iter = 1
    while iter <= maxiter && !isempty(Q)
        i,j = dequeue!(Q)
        adj1 = adjnodes(M.G1,i)
        adj2 = adjnodes(M.G2,j)
        M.f[i] = j
        push!(L1,i)
        push!(L2,j)
        # add vote to adjacent nodes of both networks
        Ddiff = wavescorevote(M.phi, i,j, adj1,adj2, adjcount1,adjcount2)
        M.D[adj1,adj2] += Ddiff

        if approxobjval # not really the obj. val. but it gives you an idea
            a1 = collect(intersect(adj1,L1))
            a2 = M.f[a1]
            aix = findin(a2,intersect(adj2,L2))
            a1 = a1[aix]
            a2 = a2[aix]
            M.objval += sum(Ddiff[sub2ind(size(Ddiff),findin(adj1,a1),
                                          findin(adj2,a2))])/min(nnz(M.G1),nnz(M.G2))
            M.objval += M.D[i,M.f[i]]/min(nnz(M.G1),nnz(M.G2))
            push!(M.history, M.objval)
        end

        updatepriorityqueue!(M,Q,L1,L2,i,j,adj1,adj2)

        print("\rIteration $iter/$maxiter")
        if approxobjval
            print(", approx. of objective value: $(M.objval)")
        end
        flush(STDOUT)
        iter += 1
    end
    if approxobjval
        println("\nPlot of iteration number vs. approx. of objective value:")
        print(lineplot(M.history))
        print("Approx. of objective value: $(M.objval)")
    end
    x = measure(M.phi,M.f)
    M.objval, M.objval_edge, M.objval_node = score(x), score(x.s), score(x.t)
    println("\nObjective value: $(M.objval), edge $(M.objval_edge), node $(M.objval_node)")
    M
end

function wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
              S::AbstractMatrix, beta::Float64,
              seeds=zeros(Int,(0,0));skipalign=false)
    M = WaveModel(G1,G2, ConvexCombMeasure(WECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    !skipalign && align!(M,approxobjval=false)
    M
end

function dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                  S::AbstractMatrix, beta::Float64,
                  seeds=zeros(Int,(0,0));skipalign=false)
    M = WaveModel(G1,G2, ConvexCombMeasure(DWECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    !skipalign && align!(M,approxobjval=false)
    M
end

# this shuffles each input network independently before passing to wave
# this is to get rid of node order biases when you have the same node set
# i.e., if both networks have the same ordered array of nodes
# then the permutation 1:n is much more likely to be the resultant
# alignment than not
# this is for use in evaluation, where we obv. do not want biases like these
function waveshuffle(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                     S::AbstractMatrix, beta::Float64,
                     method=wave,seeds=zeros(Int,(0,0)))
    p1 = randperm(size(G1,1))
    p2 = randperm(size(G2,1))
    q1 = invperm(p1)
    q2 = invperm(p2)
    H1 = G1[p1,p1]
    H2 = G2[p2,p2]

    I,J = ind2sub(size(Snode), 1:length(Snode))
    T = reshape(S[sub2ind(size(S), p1[I], p2[J])], size(S))

    M = method(H1,H2,T,beta,seeds)
    N = method(G1,G2,S,beta,seeds,skipalign=true)
    N.f = p2[M.f][q1]
    x = measure(N.phi,N.f)
    N.objval,N.objval_edge,N.objval_node = score(x),score(x.s),score(x.t)
    N
end
