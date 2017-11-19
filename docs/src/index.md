# Introduction

[DynaWAVE](https://github.com/vvjn/DynaWAVE.jl) is our software tool
for pairwise global alignment of dynamic
networks. Dynamic networks are networks that evolve over time. Until
recently, existing methods for network alignment (NA) were limited to being
able to only align static networks. However, most complex real-world
systems evolve over time and should thus be modeled as dynamic
networks. Thus, [DynaMAGNA++](https://www3.nd.edu/~cone/DynaMAGNA++/)
[1] was introduced as a proof-of-concept method for aligning dynamic
networks. However, DynaMAGNA++ does not necessarily scale well to
larger networks.

DynaWAVE [2], on the other hand, scales well to larger networks in terms
of alignment quality and runtime compared to DynaMAGNA++. While
DynaWAVE is less accurate but faster than DynaMAGNA++ for smaller
networks, DynaWAVE is more accurate and faster than DynaMAGNA++ for
larger networks.
DynaWAVE is an extension of WAVE [3]. While both WAVE and DynaWAVE
optimize edge as well as node conservation across the aligned
networks, WAVE conserves static edges and similarity between static
node neighborhoods, while DynaWAVE conserves dynamic edges (events)
and similarity between evolving node neighborhoods.

An application domain where DynaWAVE and DynaMAGNA++ are useful is
computational biology. DynaWAVE can be used for alignment of protein
interaction networks that evolve over time [4].

The following sections contain a tutorial on how to align dynamic
networks using DynaWAVE, and how to align static networks using WAVE.

[1] Alignment of dynamic networks, V. Vijayan, D. Critchlow, and
T. Milenković, Bioinformatics, Volume 33, Issue 14, 15 July 2017,
Pages i180–i189, (<https://doi.org/10.1093/bioinformatics/btx246>).

[2] Aligning dynamic networks with DynaWAVE, V. Vijayan and
T. Milenković (2017), under revision.

[3] Simultaneous Optimization of Both Node and Edge Conservation in
Network Alignment via WAVE, Y. Sun, J. Crawford, J. Tang, and
T. Milenković, in Proceedings of the Workshop on Algorithms in
Bioinformatics (WABI), Atlanta, GA, USA, September 10-12, 2015, pages
16-39 (<https://doi.org/10.1007/978-3-662-48221-6_2>).

[4] The post-genomic era of biological network alignment, EURASIP
Journal on Bioinformatics and Systems Biology, December 2015, 2015:3,
F. E. Faisal, L. Meng, J. Crawford, and T. Milenković (<https://doi.org/10.1186/s13637-015-0022-9>).

# Installation

DynaWAVE requires [Julia 0.6](https://julialang.org/).
DynaWAVE can be installed by going to the base directory of this
software and running `julia install_dynawave.jl` in the command-line.
After this, we can start `julia` and use the DynaWAVE module.
To follow the instructions below more easily, start julia in the
`examples/` directory.

We load the `DynaWAVE` module to align networks, as well as the `NetalignUtils` and
`NetalignMeasures` modules to read network files and run  related functions.

```julia
using DynaWAVE, NetalignUtils, NetalignMeasures
```

# Aligning two dynamic networks

## A toy example

First, we read in the dynamic networks `ev1.txt` and `ev2.txt` (from the
`examples/` directory). Note that the two networks are very similar,
with corresponding nodes (e.g. `g1a` corresponds to `g2a`, and so on).

The dynamic networks are represented in the events list format. That
is, each line represents and event, which is a interaction between two
nodes from the start time to the end time. For example, on the first
line of `ev1.txt`, the event represents an interaction between nodes `g1a`
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
R = readlistmat("evsim.txt", net1.nodes, net2.nodes)
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
f = dynawave(net1.G, net2.G, R)
```

`f` contains the alignment between the two networks. We construct the aligned node pairs as follows.

```julia
nodepairs = hcat(net1.nodes, net2.nodes[f])
```

We can write the alignment to file as follows.

```julia
writedlm("exalnfile.txt", nodepairs)

```

## Computing node similarities

### Topological node similarities

Dynamic GDVs [5], i.e.  dynamic graphlet degree vectors, are node
descriptors that take both the local structural topology and local
temporal topology of nodes in a dynamic network into account. The dynamic GDVs
of a network can be calculated using the dynamic graphlet counting
code available [here](https://www3.nd.edu/~cone/DG/). In the following
example, we will align two networks from the `examples/` directory, a
dynamic yeast network, and the same network with 10% of its events
randomized, while using dynamic GDVs as node similarities.

First, we read the networks and the dynamic GDVs into memory.

``` julia
net1 = readeventlist("yeastlc_original_tw_1_1.dy")
dgdv1 = readgdv("yeastlc_original_dgdv_6_4_1.txt", net1.nodes)
net2 = readeventlist("yeastlc_rnd_0.10_1_tw_1_1.dy")
dgdv2 = readgdv("yeastlc_rnd_0.10_1_dgdv_6_4_1.txt", net2.nodes)
```

Second, we calculate node similarities between node pairs in the two
networks using dynamic GDVs, making use of the PCA-based technique
described in [5].

``` julia
S = NodeSimMeasure(:pcagdvs,dgdv1,dgdv2).S
```

Finally, we align the two networks using DynaWAVE.

``` julia
f = dynawave(net1.G, net2.G, S)
```

We can construct the node pairs and calculate node correctness as follows

``` julia
nodepairs = hcat(net1.nodes, net2.nodes[f])

nc = mean(net1.nodes .== net2.nodes[f])
```

[5] Exploring the structure and function of temporal networks with
dynamic graphlets, Y. Hulovatyy,  H. Chen, and  T. Milenković,
Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages i171–i180,
(<https://doi.org/10.1093/bioinformatics/btv227>).

### External node similarities

Instead of dynamic GDVs, we can alternatively use BLAST E-values to
align networks as follows.

First, we read in the E-values to memory. We will be using the same
yeast networks as in the above sub-section.

``` julia
E = readevalues("yeastlc_yeastlc_evalues.txt", net1.nodes, net2.nodes)
```

Then, we convert the E-values to node similarities. This is done by
converting each E-value to -log(E-value), and then dividing the
resulting values with the maximum.

``` julia
S = NodeSimMeasure(:evalues, E).S
```

Finally, we align the two networks.

``` julia
f = dynawave(net1.G, net2.G, S)
```

## Creating network instances from network models

We can create random network instances from three network models as in
the DynaWAVE paper. We create a random instance of a 1000-node dynamic
network using the GEO-GD network model, with parameter $p = 0.3$, and
linear node arrival. We let the timespan range from 0 to 30 seconds,
initializing the network with a 5-node clique. 

```julia
G = rand(GEOGD(0.3, 1, :linear), 1000, 30, 5)
```

Here, we create a random instance of a 1000-node dynamic
network using the SF-GD network model, with parameters $p = 0.3, q = 0.7$, and
exponential node arrival. We let the timespan range from 0 to 30 seconds,
initializing the network with a 5-node clique. 

```julia
G = rand(SFGD(0.3, 0.7, :exp), 1000, 30, 5)
```

Here, we create a random instance of a 1000-node dynamic
network using the SNE network model, with parameters $\lambda =
0.032, \alpha = 0.8, \beta = 0.002$, and
quadratic node arrival. We let the timespan range from 0 to 30 seconds,
initializing the network with a 5-node clique. 

```julia
G = rand(SocialNE(0.032, 0.8, 0.002, :quad), 1000, 30, 5)
```

## Adding noise to a dynamic network

Given a dynamic network, here we add noise to the network. There are
two methods we use to add noise to a dynamic network: a strict version
(`strict_events_shuffle`) that only changes the event times in the network (from page 15/30 of
[6]), and a non-strict version (`links_shuffle`) that change both event times and links
between nodes in the network (from page 16/39 of [6]).

Here, we add 30% noise to network `G` using the strict version.

```julia
net1 = readeventlist("yeastlc_original_tw_1_1.dy")
G = net1.G

G30 = strict_events_shuffle(G, 0.30)
```

Here, we add 30% noise to network `G` using the non-strict version.

```julia
net1 = readeventlist("yeastlc_original_tw_1_1.dy")
G = net1.G

G30 = links_shuffle(G, 0.30)
```

[6] Modern temporal network theory: a colloquium, Petter Holme,
European Physical Journal B (2015),
(<https://doi.org/10.1140/epjb/e2015-60657-4>).

## Converting a dynamic network to a static network

We can flatten a dynamic network to a static network as follows. There
is an edge between two nodes in the static network if there is atleast
one interaction (dynamic edge) between two nodes in the dynamic
network.

```julia
net1 = readeventlist("yeastlc_original_tw_1_1.dy")
G = net1.G

Gstatic = flatten(G)
```

# Aligning two static networks

## A toy example

This software tool can also be used to align two static networks using
the WAVE method [3]. First, we read
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
f = wave(net1.G, net2.G, R)
```

Same as before, `f` contains the alignment between the two networks and we construct the aligned node pairs as follows.

```julia
nodepairs = hcat(net1.nodes, net2.nodes[f])
```

We can write the alignment to file as follows.

```julia
writedlm("exalnfile.txt", nodepairs)
```

## Computing node similarities

### Topological node similarities

Static GDVs, or GDVs, are node descriptors that take the local
topology of nodes in a static network into account. The GDVs of a
network can be calculated using
[GraphCrunch](https://en.wikipedia.org/wiki/GraphCrunch). In the
following example, we will align, from the `examples/` directory, a
yeast network to itself.

First, we read the networks and the GDVs into memory.

``` julia
net1 = readgw("0Krogan_2007_high.gw")
gdv1 = readgdv("0Krogan_2007_high.ncount.ndump2", net1.nodes)
net2 = readgw("0Krogan_2007_high.gw")
gdv2 = readgdv("0Krogan_2007_high.ncount.ndump2", net2.nodes)
```

Second, we calculate node similarities between node pairs in the two
networks using GDV similarity.

``` julia
S = NodeSimMeasure(:gdvs, gdv1, gdv2).S
```

Finally, we align the two networks using WAVE.

``` julia
f = wave(net1.G, net2.G, S)
```

We can construct the node pairs and calculate node correctness as follows

``` julia
nodepairs = hcat(net1.nodes, net2.nodes[f])

nc = mean(net1.nodes .== net2.nodes[f])
```

### External node similarities

Of course, we can similarly use E-values to align two networks as
well. We will be using the same yeast networks as in the above
sub-section. We convert the E-values to node similarities and then align.

``` julia
E = readevalues("yeastlc_yeastlc_evalues.txt", net1.nodes, net2.nodes)
S = NodeSimMeasure(:evalues, E).S
f = wave(net1.G, net2.G, S)
```

# Comparison of static NA (via WAVE) and dynamic NA (via DynaWAVE)

In the DynaWAVE paper, we compared DynaWAVE to WAVE (and MAGNA++ and
DynaMAGNA++). The following is an example demonstrating how to make
this comparison.

First, we load the dynamic networks, their corresponding dynamic GDVs,
and calculate node similarities.

``` julia
t1 = readeventlist("yeastlc_original_tw_1_1.dy")
dgdv1 = readgdv("yeastlc_original_dgdv_6_4_1.txt", t1.nodes)

t2 = readeventlist("yeastlc_rnd_0.10_1_tw_1_1.dy")
dgdv2 = readgdv("yeastlc_rnd_0.10_1_dgdv_6_4_1.txt", t2.nodes)

Sdgdv = NodeSimMeasure(:pcagdvs,dgdv1,dgdv2).S
```

Then, we load the corresponding flattened static networks, their
corresponding GDVs, and calculate node similarities calculated using
the PCA-based technique described in [5].

``` julia
s1 = readgw("yeastlc_original.gw")
gdv1 = readgdv("yeastlc_original.gdv.ndump2", s1.nodes)

s2 = readgw("yeastlc_original.gw")
gdv2 = readgdv("yeastlc_original.gdv.ndump2", s2.nodes)

Sgdv = NodeSimMeasure(:pcagdvs,gdv1,gdv2).S
```

Notice that the flattened version of both `t1` and `t1` result in the
same network as `s1` (modulo node permutations) due to our randomization model. That is:

```julia
flatten(t1[sortperm(t1.nodes)]) == flatten(t2[sortperm(t2.nodes)]) == s1[sortperm(s1.nodes)]
# true
```

Then, we align the networks using both WAVE and DynaWAVE. Here, we set
the $\beta$ parameter to `0.5` as in the paper.

``` julia
fdynamic = dynawave(t1.G, t2.G, Sdgdv, 0.5)

fstatic = wave(s1.G, s2.G, Sgdv, 0.5)
```

We can calculate node correctness as follows.

``` julia
nc_dynamic = nodecorrectness(fdynamic, t1.nodes, t2.nodes)

nc_static = nodecorrectness(fstatic, s1.nodes, s2.nodes)
```

Finally, above, we aligned the yeast dynamic network to the same network with
10% randomization. We can also align the yeast dynamic network to itself.

``` julia
S = NodeSimMeasure(:pcagdvs,gdv1,gdv1).S
f = dynawave(t1.G, t1.G, S, 0.5)
nodecorrectness(f, t1.nodes, t1.nodes)
```
