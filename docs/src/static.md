# Static NA

# [Running WAVE with pre-computed node similarities](@id sec6)

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

# Computing node similarities and then running WAVE

## [Computing topological node similarities](@id sec7)

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

## [Computing external node similarities from BLAST E-values](@id sec8)

Of course, we can similarly use E-values to align two networks as
well. We will be using the same yeast networks as in the above
sub-section. We convert the E-values to node similarities and then align.

``` julia
E = readevalues("yeastlc_yeastlc_evalues.txt", net1.nodes, net2.nodes)
S = NodeSimMeasure(:evalues, E).S
f = wave(net1.G, net2.G, S)
```

# [Comparison of static and dynamic NA](@id sec9)

In the DynaWAVE paper, we compared DynaWAVE to WAVE (and MAGNA++ and
DynaMAGNA++). The following is an example demonstrating how to make
this comparison using WAVE and DynaWAVE.

First, we load the two dynamic networks, their corresponding dynamic GDVs,
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

We can flatten a dynamic network to a static network as follows. There
is an edge between two nodes in the static network if there is atleast
one interaction (dynamic edge) between two nodes in the dynamic
network.

```julia
net1 = readeventlist("yeastlc_original_tw_1_1.dy")
Gdynamic = net1.G

Gstatic = flatten(Gdynamic)
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
