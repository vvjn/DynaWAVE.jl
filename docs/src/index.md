# Tutorial

# I. Introduction

DynaWAVE is our software tool for pairwise global alignment of dynamic networks.

The following sections contain a tutorial on how to align dynamic
networks using DynaWAVE.

Namely, Section [II](@ref secm2) contains installation
instructions. Section [III](@ref secm3) describes how to align two
dynamic networks when pre-computed node similarities are available as
input to DynaWAVE. On the other hand, Section [IV](@ref secm4)
describes how to align two dynamic networks without having available
pre-computed node similarities. Section [V](@ref secm5) describes how
to mimic some of the experiments from our DynaWAVE paper.

[1] Aligning dynamic networks with DynaWAVE, V. Vijayan and
T. Milenković (2017), Bioinformatics (2017), btx841,
(<https://doi.org/10.1093/bioinformatics/btx841>).

# [II. Installation](@id secm2)

[Julia 0.6](https://julialang.org/) needs to be installed prior to
using DynaWAVE. DynaWAVE can be installed by starting Julia and executing the command `Pkg.add("DynaWAVE")`.


# [III. Running DynaWAVE with pre-computed node similarities](@id secm3)

If you have pre-computed node similarities on which you wish to run
DynaWAVE, then you may perform the following instructions. If you wish to
compute node similarities as in the DynaWAVE paper, then go to Section
[IV](@ref secm4).

**1.** Start by navigating to the `bin/` directory from the base
directory of this package (`Pkg.dir("DynaWAVE")`), and running the command `julia
dynawave.jl`. This shows examples of how to align two dynamic
networks, given node similarities between them. Below are detailed
instructions on how to align two example networks.

**2.** Choose two networks to align. DynaWAVE accepts networks in the
event list format (see an [example network 1](ev1.txt) and an [example
network 2](ev2.txt)). Each line in the event list format contains an
event, which consists of the event's start time, end time, first node,
and second node, respectively.

**3.** Choose a node similarity file. DynaWAVE accepts node
similarities in the node pair list format (see an [example node
similarity file](evsim.txt) corresponding to the above two
networks). The node similarity file consists of three columns, the
first column being nodes from network 1, the second column being nodes
from network 2, and the third column being the node similarities. Each
node similarity must be between 0 and 1, inclusive. The node
similarity file does not necessarily need to contain similarities
between all pairs of nodes across the networks.

**4.** Choose the output alignment file name. DynaWAVE will output the
alignment file to this file.

**5.** Given network 1 file, `ev1.txt`, network 2 file, `ev2.txt`,
node similarity file, `evsim.txt`, and the output alignment file
`output_alignment.txt`, the following command is run in order to align
the two networks:

```
julia dynawave.jl ev1.txt ev2.txt evsim.txt output_alignment.txt
```

The above command will produce the file `output_alignment.txt`, which
contains the alignment. It will also print various statistics related
to the alignment.


# [IV. Computing node similarities and then running DynaWAVE](@id secm4)

To follow the instructions below, start `julia` in the `examples/`
directory from the base directory of this software, and run the
following command to load the DynaWAVE library/module:

```julia
using DynaWAVE, NetalignUtils, NetalignMeasures
```

If you have pre-computed node similarities on which you wish to run
DynaWAVE, then you can use the interface at Section [III](@ref
secm3). If you wish to compute the topological node similarities as in
the DynaWAVE paper, then perform the instructions in the following
[sub-section](@ref sec2). If you wish to compute node similarities
from BLAST E-values to align dynamic protein interaction networks,
then go to the [sub-section](@ref sec3) after that.

### [Computing topological node similarities](@id sec2)

Dynamic GDVs [2], i.e., dynamic graphlet degree vectors, are node
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
described in [2].

``` julia
S = NodeSimMeasure(:pcagdvs,dgdv1,dgdv2).S
```

(Optionally, you can run `writelistmat("nodesims.txt", S, net1.nodes,
net2.nodes)` to save the node similarities to a file, and you can run `S
= readlistmat("nodesims.txt", net1.nodes, net2.nodes)`, to read the
node similarities from the file.)

Finally, we align the two networks using DynaWAVE.

``` julia
f = dynawave(net1.G, net2.G, S)
```

`f` contains the alignment between the two networks. We can construct the node pairs and calculate node correctness as follows

``` julia
nodepairs = hcat(net1.nodes, net2.nodes[f])

nc = mean(net1.nodes .== net2.nodes[f])
```

We can write the alignment to file as follows.

```julia
writedlm("yeastlc_aln.txt", nodepairs)

```

[2] Exploring the structure and function of temporal networks with
dynamic graphlets, Y. Hulovatyy,  H. Chen, and  T. Milenković,
Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages i171–i180,
(<https://doi.org/10.1093/bioinformatics/btv227>).

### [Computing external node similarities from BLAST E-values](@id sec3)

Instead of dynamic GDVs, we can alternatively use BLAST E-values to
align protein interaction networks as follows.

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

# [V. Mimicking experiments from DynaWAVE paper](@id secm5)

If you wish to create synthetic dynamic networks from network models
as in the DynaWAVE paper, then go the following [sub-section](@ref
sec4). If you wish to add noise to a dynamic network as in the
DynaWAVE paper, then go to the [sub-section](@ref sec5) after that.

## [Creating synthetic dynamic networks](@id sec4)

In this sub-section, we create random network instances from three
network models as in the DynaWAVE paper. Here, we create a random
instance of a 1000-node dynamic network using the GEO-GD network
model, with parameter $p = 0.3$, and linear node arrival. We let the
timespan range from 0 to 30 seconds, initializing the network with a
5-node clique.

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

## [Adding noise to a dynamic network](@id sec5)

In this sub-section, given a dynamic network, we add noise to the
network. There are two methods we use in the DynaWAVE paper to add
noise to a dynamic network: a strict version (`strict_events_shuffle`)
that only changes the event times in the network (from page 15/30 of
[3]), and a non-strict version (`links_shuffle`) that change both
event times and links between nodes in the network (from page 16/39 of
[3]).

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

[3] Modern temporal network theory: a colloquium, Petter Holme,
European Physical Journal B (2015),
(<https://doi.org/10.1140/epjb/e2015-60657-4>).


