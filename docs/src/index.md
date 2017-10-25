# Introduction

DynaWAVE is our software tool for pairwise global alignment of dynamic
networks. Dynamic networks are networks that evolve over time. Until
recently, existing methods for network alignment were limited to being
able to only align static networks. However, most complex real-world
systems evolve over time and should thus be modeled as dynamic
networks. Thus, [DynaMAGNA++](https://www3.nd.edu/~cone/DynaMAGNA++/)
[1] was introduced as a proof-of-concept
method for aligning dynamic networks. However, DynaMAGNA++ does not
necessarily scale well to larger networks.

DynaWAVE, on the other hand, scales well to larger networks in terms
of alignment quality and runtime compared to DynaMAGNA++. While
DynaWAVE is less accurate but faster than DynaMAGNA++ for smaller
networks, DynaWAVE is more accurate and faster than DynaMAGNA++ for
larger networks. An example of an application domain where DynaWAVE
and DynaMAGNA++ are useful is computational biology - DynaWAVE can be
used for alignment of protein interaction networks that evolve over
time.

DynaWAVE is an extension of WAVE. While both WAVE and DynaWAVE
optimize edge as well as node conservation across the aligned
networks, WAVE conserves static edges and similarity between static
node neighborhoods, while DynaWAVE conserves dynamic edges (events)
and similarity between evolving node neighborhoods. WAVE appears in
the following publication [2].

[1] Alignment of dynamic networks,
V. Vijayan,  D. Critchlow, and  T. Milenković,
Bioinformatics, Volume 33, Issue 14, 15 July 2017, Pages i180–i189,
(https://doi.org/10.1093/bioinformatics/btx246).

[2] Yihan Sun, Joseph Crawford, Jie Tang, and
Tijana Milenković, Simultaneous Optimization of Both Node and Edge
Conservation in Network Alignment via WAVE, in Proceedings of the
Workshop on Algorithms in Bioinformatics (WABI), Atlanta, GA, USA,
September 10-12, 2015, pages 16-39.

# Installation

DynaWAVE can be installed by starting julia and executing the
following commands. This will install DynaWAVE and packages that it depends on.

```julia
Pkg.clone("https://github.com/vvjn/NetalignMeasures.jl")
Pkg.clone("https://github.com/vvjn/NetalignUtils.jl")
Pkg.clone("https://github.com/vvjn/DynaWAVE.jl")
```

# Aligning two dynamic networks

First, we load the DynaWAVE module as well as the NetalignUtils to
read network files.

```julia
using DynaWAVE
using NetalignUtils
```

Then, we read in the dynamic networks `ev1.txt` and `ev2.txt` (from the
`examples/` directory). Note that the two networks are very similar,
with corresponding nodes (e.g. `g1a` corresponds to `g2a`, and so on).

The dynamic networks are represented in the events list format. That
is, each line represents and event, which is a interaction between two
nodes from the start time to the end time. For example, on the first
line of `ev1.txt`e, the event represents an interaction between nodes `g1a`
and `g1b` from start time 3 to end time 5.

```julia
net1 = readeventlist("ev1.txt")
net2 = readeventlist("ev2.txt")
```
```
shell> cat ev1.txt
3 5 g1a g1b
6 7 g1a g1b
0 2 g1a g1b
0 4 g1b g1c
5 8 g1b g1c
0 8 g1c g1d
1 3 g1c g1e
5 8 g1c g1e

shell> cat ev2.txt
0 2 g2a g2b
3 5 g2a g2b
6 7 g2a g2b
0 4 g2b g2c
5 8 g2b g2c
0 8 g2c g2d
1 3 g2c g2e
5 8 g2c g2e
0 8 g2e g2f
```

Then, we read in the node similarities between the two networks from evsim.txt. Node similarities
can be generated in many ways including BLAST E-values, (dynamic) graphlet degree vector similarities, etc.

```julia
R = readmat("evsim.txt", net1.nodes, net2.nodes)
```
```
shell> cat evsim.txt
g1c g2a 0.7028744190938785
g1c g2b 0.9503056260723859
g1c g2f 0.7028744190938785
g1c g2e 0.9503056260723859
g1c g2c 1.0
g1c g2d 0.5051375860429862
g1b g2a 0.5458846833524831
[...]
g1d g2c 0.4117502719875372
g1d g2d 0.24794887906050478
```

Finally, we align the two networks using the node similarities.

```julia
f = dynawave(net1.G, net2.G, R, 0.5);
```

`f` contains the alignment between the two networks. We construct the aligned node pairs as follows.

```julia
nodepairs = hcat(net1.nodes, net2.nodes[f])
```

We can write the alignment to file as follows.

```julia
writedlm("exalnfile.txt", nodepairs)
```

# Aligning two static networks

DynaWAVE can also be used to align two static networks. First, we read
in the static networks `ex1.gw` and `ex2.gw` (from the `examples/`
directory). The two networks are very similar, with corresponding
nodes of the same names. The networks are in the LEDA format (hence,
the `.gw` extension) and so we use the `readgw` function to load it.

```julia
net1 = readgw("ex1.gw")
net2 = readgw("ex2.gw")
```
```
shell> cat ex1.gw
LEDA.GRAPH
void
void
-2
5
|{A}|
|{B}|
|{C}|
|{D}|
|{E}|
5
1 2 0 |{}|
2 3 0 |{}|
3 4 0 |{}|
4 1 0 |{}|
4 5 0 |{}|

shell> cat ex2.gw
LEDA.GRAPH
void
void
-2
6
|{A}|
|{B}|
|{C}|
|{D}|
|{E}|
|{F}|
6
1 2 0 |{}|
2 3 0 |{}|
3 4 0 |{}|
4 1 0 |{}|
4 5 0 |{}|
5 6 0 |{}|
```

Then, we read in the node similarities between the two networks from
`evgdvsim.txt`. The node similarities are also stored in a different
format and we use the `readdlm` function to load it.

```julia
R = readdlm("exgdvsim.txt", header=true)[1]
```
```
shell> cat exgdvsim.txt
5 6
0.972933 0.93658 0.972933 0.925753 0.935525 0.903174
0.939805 0.988206 0.939805 0.91356 0.933523 0.937895
0.972933 0.93658 0.972933 0.925753 0.935525 0.903174
0.917616 0.925602 0.917616 0.973495 0.910664 0.896031
0.922085 0.933688 0.922085 0.879783 0.951919 0.954437
```

Finally, we align the two networks using the node similarities.

```julia
f = wave(net1.G, net2.G, R, 0.5);
```

Same as before, `f` contains the alignment between the two networks and we construct the aligned node pairs as follows.

```julia
nodepairs = hcat(net1.nodes, net2.nodes[f])
```

We can write the alignment to file as follows.

```julia
writedlm("exalnfile.txt", nodepairs)
```
