using DataStructures
using UnicodePlots

export wave, dynawave, shufflealign, WaveModel, align!

type WaveModel{E,T<:AbstractMatrix,M<:NetalignMeasure}
    G1 :: SparseMatrixCSC{E,Int}
    G2 :: SparseMatrixCSC{E,Int}
    D :: T # vote matrix
    f :: Vector{Int} # mapping
    phi :: M # objective value
    seeds :: Vector{Tuple{Int,Int}} # initial seeds k x 2 int array
    objval :: Float64 # objective value
    es
    history :: Vector{Float64} # approx. obj. val. at each iteration

    function WaveModel{E,T,M}(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                       D::AbstractMatrix,
                       meas::NetalignMeasure,
                       seeds) where {E,T,M}
        new(G1,G2, D,
            zeros(Int,size(G1,1)),
            meas,seeds,0.0,nothing,Float64[])
    end
end
"""
    WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
              meas::NetalignMeasure, seeds = [])

Creates a problem object given two networks to align by
optimizing a `NetalignMeasure`. See [`align!`](@ref) on how perform
the alignment.

# Arguments
- `G1`,`G2` : input networks in sparse matrix format
- `meas` : `NetalignMeasure` object (see
  [NetalignMeasures.jl](https://github.com/vvjn/NetalignMeasures.jl)). WaveModel
  can optimize the DWEC measure (`DWECMeasure`), a dynamic network
  alignment measure and the WEC measure (`WECMeasure`), a static
  network alignment method.
"""
function WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                   meas::NetalignMeasure, seeds = Vector{Tuple{Int,Int}}())
    D = wavescorematrix(meas)
    WaveModel{eltype(G1),typeof(D),typeof(meas)}(G1,G2,D,meas,seeds)
end

function wavestep!(M::WaveModel,i::Int,j::Int,
                   L1::Set{Int},L2::Set{Int},
                   adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int})
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
    adj1,adj2,Ddiff
end

function initializewave!(M::WaveModel;verbose=false)
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    # adjcount[i] is the num. of times that a neighbor of i is aligned
    adjcount1 = zeros(Int,n1)
    adjcount2 = zeros(Int,n2)
    M.objval = 0.0
    L1 = Set{Int}()
    L2 = Set{Int}() # nodes already aligned
    for k = 1:length(M.seeds)
        i,j = M.seeds[k]
        adj1,adj2,Ddiff = wavestep!(M,i,j,L1,L2,adjcount1,adjcount2)

        if verbose
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
    Q
end

"""
    align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-length(M.seeds)),
                verbose=false) -> f

Performs alignment with respect to the parameters in `WaveModel`.

# Arguments
- `maxiter` : Stop aligning after `maxiter` iterations. Each iteration
    creates one aligned node pair.
- `verbose` : Prints debugging outputs
"""
function align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-length(M.seeds)),
                verbose=false)
    # intialize alignment using the seeds aligned pairs
    L1,L2,adjcount1,adjcount2 = initializewave!(M,verbose=verbose)

    # initialize priority queue
    println("Building priority queue")
    Q = buildpriorityqueue(M,L1,L2)

    iter = 1
    while iter <= maxiter && !isempty(Q)
        i,j = dequeue!(Q)
        adj1,adj2,Ddiff = wavestep!(M,i,j,L1,L2,adjcount1,adjcount2)

        if verbose # not really the obj. val. but it gives you an idea
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
        if verbose
            print(", approx. of objective value: $(M.objval)")
        end
        iter += 1
    end
    if verbose
        println("\nPlot of iteration number vs. approx. of objective value:")
        print(lineplot(M.history))
        print("Approx. of objective value: $(M.objval)")
    end
    x = measure(M.phi,M.f)
    M.objval = score(x)
    M.es = x
    println("\nObjective value: $(M.objval)")
    f
end

"""
    wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
         S::AbstractMatrix, [ beta=0.5,seeds=[] ]);
         skipalign=false,details=false) -> f [, M]

Given two static networks and node similarities between them, align the two
networks by running the WAVE algorithm (Yihan Sun, Joseph Crawford,
Jie Tang, and Tijana Milenkovic, Simultaneous Optimization of Both
Node and Edge Conservation in Network Alignment via WAVE, in
Proceedings of the Workshop on Algorithms in Bioinformatics (WABI),
Atlanta, GA, USA, September 10-12, 2015, pages 16-39).

# Arguments
- `G1`,`G2` : input networks in sparse matrix format
- `S` : node similarities beween the two networks
- `beta` : weighs between edge and node conservation, `beta=1.0` weighs
    edge conservation highly while `beta=0.0` weighs node conservation highly
- `seeds` : seed aligned node pairs. For example, if we know that the
    3rd node in the first network is aligned to the 7th node in the
    second network, and the 5th node in the first network is aligned to the
    9th node in the second network, then we will set `seeds = [(3,7), (5,9)]`

# Keyword arguments
- `skipalign` : Don't align; just return the alignment details in [`WaveModel`](@ref).
- `details` : Returns `M`, the [`WaveModel`](@ref).
"""
function wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
              S::AbstractMatrix, beta=0.5,
              seeds=Vector{Tuple{Int,Int}}();skipalign=false,details=false)
    M = WaveModel(G1,G2, ConvexCombMeasure(WECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    if !skipalign
        align!(M,verbose=false,details=details)
    else
        M.f,M
    end
end

"""
    dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
         S::AbstractMatrix, [ beta=0.5,seeds=[] ]);
         skipalign=false,details=false) -> f [, M]

Given two dynamic networks and node similarities between them, align the two
networks by running the DynaWAVE algorithm.

# Arguments
- `G1`,`G2` : input networks in sparse matrix format
- `S` : node similarities beween the two networks
- `beta` : weight between edge and node conservation, `beta=1.0` weighs
    edge conservation highly while `beta=0.0` weighs node conservation highly
- `seeds` : seed aligned node pairs. For example, if we know that the
    3rd node in the first network is aligned to the 7th node in the
    second network, and the 5th node in the first network is aligned to the
    9th node in the second network, then we will set `seeds = [(3,7), (5,9)]`

# Keyword arguments
- `skipalign` : Don't align; just return the alignment details in [`WaveModel`](@ref).
- `details` : Returns `M`, the [`WaveModel`](@ref).
"""
function dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                  S::AbstractMatrix, beta=0.5,
                  seeds=Vector{Tuple{Int,Int}}();skipalign=false,details=false)
    M = WaveModel(G1,G2, ConvexCombMeasure(DWECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    if !skipalign
        align!(M,verbose=false,details=details)
    else
        M.f,M
    end
end

"""
    shufflealign(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                      S::AbstractMatrix, beta::Float64,
                      seeds=[], details=false,method=dynawave) -> f [, M]

Given two networks and node similarities between them, before aligning
using either [`dynawave`](@ref) and [`wave`](@ref), this shuffles each
input network independently before passing the networks to `dynawave`
or `wave`. This is to get rid of node order biases when you have the
same node set; i.e., if both networks have the same ordered array of
nodes then the permutation 1:n is much more likely to be the resultant
alignment than not. This is for use in evaluation, where biases like
this tend to show much better results that is actually possible with a
particular method. We obviously do not want biases like these during evaluation.
For example, when aligning a network to itself evaluate a network alignment
method, the node order of the two networks will be the same. This
results in the network alignment method "knowing" the true node mapping,
and often producing alignments of much higher quality that is actually
possible with the method were the true node mapping not known.    

# Arguments
- `method` : If `method = dynawave`, run the DynaWAVE algorithm after
    shuffling the networks. If `method = wave`, run the WAVE algorithm
    after shuffling the networks.
- See [`dynawave`](@ref) and [`wave`](@ref) for the other arguments.
    
# Output
- `f` : The resulting alignment has the correct node order with respect to the input `G1` and `G2` networks; i.e., it is "unshuffled".    
"""    
function shufflealign(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                      S::AbstractMatrix, beta::Float64,
                      seeds=Vector{Tuple{Int,Int}}();
                      details=false,method=dynawave)
    p1 = randperm(size(G1,1))
    p2 = randperm(size(G2,1))
    q1 = invperm(p1)
    q2 = invperm(p2)
    H1 = permute(G1,p1,p1)
    H2 = permute(G2,p2,p2)

    I,J = ind2sub(size(S), 1:length(S))
    T = reshape(S[sub2ind(size(S), p1[I], p2[J])], size(S))

    _,M = method(H1,H2,T,beta,seeds,details=true)
    _,N = method(G1,G2,S,beta,seeds,skipalign=true,details=true)
    N.f = p2[M.f][q1]
    x = measure(N.phi,N.f)
    N.objval = score(x)
    N.es = x
    if details N.f,N else N.f end
end
