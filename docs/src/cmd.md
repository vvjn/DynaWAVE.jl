# Command-line tool

Here are instructions on running the command-line tool `bin/dynawave.jl`.

Start DynaWAVE by running the command `julia bin/dynawave.jl`. Examples of how to run DynaWAVE can be found by running `julia bin/dynawave.jl --help`. Below we show how to align two example networks.

Choose two networks to align. DynaWAVE accepts networks in the event list format (see an example network 1 and an example network 2). Each line in the event list format contains an event, which consists of the event's start time, end time, first node, and second node, respectively.

Choose a node similarity file. DynaWAVE accepts node similarities in the node pair list format (see an example node similarity file corresponding to the above two networks). The node similarity file consists of three columns, the first column being nodes from network 1, the second column being nodes from network 2, and the third column being the node similarities. Each node similarity must be between 0 and 1, inclusive. The node similarity file does not necessarily need to contain similarities between all pairs of nodes across the networks.

Choose the output alignment file name. DynaWAVE will output the alignment file to this file.

The above three steps allows us to perform dynamic network alignment using DynaWAVE. Namely, given network 1 file, `ev1.txt`, network 2 file, `ev2.txt`, node similarity file, `exev1ev2.txt`, and the output alignment file `output_alignment.txt`, the following command is run in order to align the two networks:

```julia bin/dynawave.jl ev1.txt ev2.txt exev1ev2.txt output_alignment.txt```

The above command will produce the file `output_alignment.txt`, which contains the alignment. It will also print various statistics related to the alignment.
