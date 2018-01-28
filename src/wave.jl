export wave, dynawave, shufflealign, WaveModel, align!

type WaveModel{E,T<:AbstractMatrix,M<:NetalignMeasure}
    G1 :: SparseMatrixCSC{E,Int}
    G2 :: SparseMatrixCSC{E,Int}
    D :: T # vote matrix
    f :: Vector{Int} # mapping
    phi :: M # objective value
    seeds :: Vector{Tuple{Int,Int}} # initial seeds k x 2 int array
    objval :: Float64 # objective value
    nas # :: NetalignScore
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

Internal function that creates a problem structure given two networks
(static or dynamic) to align by and `NetalignMeasure` to optimize
over. See [`align!`](@ref) on how perform the alignment.

# Arguments
- `G1`,`G2` : input networks in sparse matrix format
- `meas` : `NetalignMeasure` object (see
  [NetalignMeasures.jl](https://github.com/vvjn/NetalignMeasures.jl)). `WaveModel`
  can optimize the DWEC measure (`DWECMeasure`), a dynamic network
  alignment measure and the WEC measure (`WECMeasure`), a static
  network alignment method.
"""
function WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                   meas::NetalignMeasure, seeds = Vector{Tuple{Int,Int}}())
    D = wavescorematrix(meas)
    WaveModel{eltype(G1),typeof(D),typeof(meas)}(G1,G2,D,meas,seeds)
end

function initializewave!(M::WaveModel;verbose=false)
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    Q = buildpriorityqueue(M)
    # adjcount[i] is the num. of times that a neighbor of i is aligned
    adjcount1 = zeros(Int,n1)
    adjcount2 = zeros(Int,n2)
    M.objval = 0.0
    L1 = Set{Int}()
    L2 = Set{Int}() # nodes already aligned
    for k = 1:length(M.seeds)
        i,j = M.seeds[k]
        wavestep!(M,i,j,L1,L2,adjcount1,adjcount2,Q,verbose)
    end
    (Q,L1,L2,adjcount1,adjcount2)
end

function buildpriorityqueue(M::WaveModel)
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    if issparse(M.D)
        I,J,V = findnz(M.D)
        kv = (((i,j),-v) for (i,j,v) in zip(I,J,V))
    else
        kv = (((i,j),-M.D[i,j]) for i=1:n1, j=1:n2)
    end
    PriorityQueue{Tuple{Int,Int},Float64}(kv)
end

function updatepriorityqueue!(M::WaveModel,Q::PriorityQueue,
                              L1::Set{Int},L2::Set{Int}, i::Int,j::Int,
                              adj1::AbstractVector{Int},adj2::AbstractVector{Int})
    n1 = size(M.G1,1); n2 = size(M.G2,1)
    # remove cross node pairs containing i or j
    if issparse(M.D)
        for ip = setdiff(1:n1,L1); if haskey(Q,(ip,j)) dequeue!(Q,(ip,j)) end end
        for jp = setdiff(1:n2,L2); if haskey(Q,(i,jp)) dequeue!(Q,(i,jp)) end end
    else
        for ip = setdiff(1:n1,L1); dequeue!(Q,(ip,j)) end
        for jp = setdiff(1:n2,L2); dequeue!(Q,(i,jp)) end
    end
    # update priority queue with the new votes
    for ip in setdiff(adj1,L1), jp in setdiff(adj2,L2)
        Q[(ip,jp)] = -M.D[ip,jp]
    end
    Q
end

function wavestep!(M::WaveModel,i::Int,j::Int,
                   L1::Set{Int},L2::Set{Int},
                   adjcount1::AbstractVector{Int},adjcount2::AbstractVector{Int},
                   Q::PriorityQueue, verbose)
    if verbose # not really the obj. val. but it gives you an idea
        M.objval += M.D[i,j]
        push!(M.history, M.objval)
    end

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
    updatepriorityqueue!(M,Q,L1,L2,i,j,adj1,adj2)
end

"""
    align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-length(M.seeds)),
                verbose=false) -> f

Internal function that performs alignment with respect to the parameters in `WaveModel`.

# Arguments
- `maxiter` : Stop aligning after `maxiter` iterations. Each iteration
    creates one aligned node pair.
- `verbose` : Prints debugging outputs
"""
function align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-length(M.seeds)),
                verbose=false)
    # intialize alignment using the seeds aligned pairs
    Q,L1,L2,adjcount1,adjcount2 = initializewave!(M,verbose=verbose)

    # initialize priority queue
    println("Building priority queue")

    iter = 1
    while iter <= maxiter && !isempty(Q)
        i,j = dequeue!(Q)
        wavestep!(M,i,j,L1,L2,adjcount1,adjcount2,Q,verbose)

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
    M.nas = x
    println("\nObjective value: $(M.objval)")
    alnfillrandom!(M.f,size(M.D,2)) # fill unfilled mappings randomly
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
- `G1`,`G2` : Input networks in sparse matrix format.
- `S` : Node similarities beween the two networks.
- `beta` : Weighs between edge and node conservation, `beta=1.0` weighs
    edge conservation highly while `beta=0.0` weighs node conservation highly
- `seeds` : Seed aligned node pairs. For example, if we know that the
    3rd node in the first network is aligned to the 7th node in the
    second network, and the 5th node in the first network is aligned to the
    9th node in the second network, then we will set `seeds = [(3,7), (5,9)]`.

# Keyword arguments
- `skipalign` : Don't align; just return the alignment details in [`WaveModel`](@ref).
- `details` : Returns `M`, the [`WaveModel`](@ref).

# Output
- `f` : Alignment, i.e. node mapping from `G1` to `G2`. `f[i] describes
    node pair `nodes1[i], nodes2[f[i]]`, where `nodes1` and `nodes2` are
    the vectors of node names corresponding to `G1` and `G2` respectively.
"""
function wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
              S::AbstractMatrix, beta=0.5,
              seeds=Vector{Tuple{Int,Int}}();skipalign=false,details=false,verbose=false)
    if issparse(S)
        S = SparseMatrixLIL(S)
    end
    M = WaveModel(G1,G2, ConvexCombMeasure(WECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    if !skipalign
        align!(M,verbose=verbose)
    end
    if details M.f,M else M.f end
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

# Output
- `f` : Alignment, i.e. node mapping from `G1` to `G2`. `f[i] describes
    node pair `nodes1[i], nodes2[f[i]]`, where `nodes1` and `nodes2` are
    the vectors of node names corresponding to `G1` and `G2` respectively.
"""
function dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                  S::AbstractMatrix, beta=0.5,
                  seeds=Vector{Tuple{Int,Int}}();skipalign=false,details=false,verbose=false)
    if typeof(S) <: SparseMatrixCSC
        S = SparseMatrixLIL(S)
    end
    M = WaveModel(G1,G2, ConvexCombMeasure(DWECMeasure(G1,G2,S),
                                           NodeSimMeasure(S),beta),seeds)
    if !skipalign
        align!(M,verbose=verbose)
    end
    if details M.f,M else M.f end
end

"""
    shufflealign(method, G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                      S::AbstractMatrix, beta::Float64,
                      seeds=[], details=false) -> f [, M]

Given two networks and node similarities between them, before aligning
using either [`dynawave`](@ref) or [`wave`](@ref), this shuffles each
input network independently before passing the networks to `dynawave`
or `wave`. This is to remove node order biases when the two networks
being aligned have similar node sets with similar node orderings. This
is for use in evaluation, where biases like this tend to show much
better results that is actually possible with a particular method. We
obviously do not want biases like these during evaluation.

For example, when aligning a network to itself evaluate a network
alignment method, the node order of the two networks will be the same
and the node set will be the same. This results in the network
alignment method "knowing" the true node mapping, and often producing
alignments of much higher quality that is actually possible with the
method were the true node mapping not known. This method randomizes
(shuffles) the node order before handing it to `dynawave` or `wave`
and then returns the deshuffled alignment.

# Arguments
- `method` : If `method = dynawave`, run the DynaWAVE algorithm after
    shuffling the networks. If `method = wave`, run the WAVE algorithm
    after shuffling the networks.
- See [`dynawave`](@ref) and [`wave`](@ref) for the other arguments.

# Output
- `f` : The resulting alignment has the correct node order with respect to the input `G1` and `G2` networks; i.e., it is "deshuffled".
"""
function shufflealign(method::Function, G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                      S::AbstractMatrix, beta::Float64,
                      seeds=Vector{Tuple{Int,Int}}();
                      details=false)
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
    N.nas = x
    if details N.f,N else N.f end
end
