# DynaWAVE

[![Build Status](https://travis-ci.org/vvjn/DynaWAVE.jl.svg?branch=master)](https://travis-ci.org/vvjn/DynaWAVE.jl) [![Coverage Status](https://coveralls.io/repos/vvjn/DynaWAVE.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/vvjn/DynaWAVE.jl?branch=master) [![codecov.io](http://codecov.io/github/vvjn/DynaWAVE.jl/coverage.svg?branch=master)](http://codecov.io/github/vvjn/DynaWAVE.jl?branch=master)

Network alignment of dynamic(temporal) networks.

# Installation

DynaWAVE can be installed by starting julia and cloning this package. This will install DynaWAVE and packages
that it depends on.

```julia
Pkg.clone("https://github.com/vvjn/DynaWAVE.jl")
```

# Aligning two dynamic networks

First, we load the required modules.

```julia
using DynaWAVE
using NetalignUtils
```

Then, we read in the dynamic networks ev1.txt and ev2.txt (from the examples/ directory). Note that
the two networks are very similar, with corresponding nodes (e.g. g1a corresponds to g2a, and so on).

```julia
net1 = readeventlist("ev1.txt")
net2 = readeventlist("ev2.txt")

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
res = dynawave(net1.G, net2.G, R, 0.5);
```

`res.f` contains the alignment between the two networks. We construct the aligned node pairs as follows.

```julia
hcat(net1.nodes, net2.nodes[res.f])
```

# Aligning two static networks

DynaWAVE can also be used to align two static networks. First, we read in the static networks ex1.gw and
ex2.gw (from the examples/ directory). The two networks are very similar, with corresponding nodes of the
same names. The networks are in the LEDA format (hence, the .gw extension) and so we use the `readgw`
function to load it.

```julia
net1 = readgw("ex1.gw")
net2 = readgw("ex2.gw")

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

Then, we read in the node similarities between the two networks from evgdvsim.txt. The node similarities
are also stored in a different format and we use the `readdlm` function to load it.

```julia
R = readdlm("exgdvsim.txt", header=true)[1]

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
res = wave(net1.G, net2.G, R, 0.5);
```

Same as before, `res.f` contains the alignment between the two networks and we construct the aligned node pairs as follows.

```julia
hcat(net1.nodes, net2.nodes[res.f])
```
