var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#I.-Introduction-1",
    "page": "Tutorial",
    "title": "I. Introduction",
    "category": "section",
    "text": "DynaWAVE is our software tool for pairwise global alignment of dynamic networks.The following sections contain a tutorial on how to align dynamic networks using DynaWAVE.Namely, Section II contains installation instructions. Section III describes how to align two dynamic networks when pre-computed node similarities are available as input to DynaWAVE. On the other hand, Section IV describes how to align two dynamic networks without having available pre-computed node similarities. Section V describes how to mimic some of the experiments from our DynaWAVE paper.[1] Aligning dynamic networks with DynaWAVE, V. Vijayan and T. Milenković (2017), under revision."
},

{
    "location": "index.html#secm2-1",
    "page": "Tutorial",
    "title": "II. Installation",
    "category": "section",
    "text": "Julia 0.6 needs to be installed prior to using DynaWAVE. DynaWAVE can be installed by starting Julia and executing the command Pkg.add(\"DynaWAVE\")."
},

{
    "location": "index.html#secm3-1",
    "page": "Tutorial",
    "title": "III. Running DynaWAVE with pre-computed node similarities",
    "category": "section",
    "text": "If you have pre-computed node similarities on which you wish to run DynaWAVE, then you may perform the following instructions. If you wish to compute node similarities as in the DynaWAVE paper, then go to Section IV.1. Start by navigating to the bin/ directory from the base directory of this package (Pkg.dir(\"DynaWAVE\")), and running the command julia dynawave.jl. This shows examples of how to align two dynamic networks, given node similarities between them. Below are detailed instructions on how to align two example networks.2. Choose two networks to align. DynaWAVE accepts networks in the event list format (see an example network 1 and an example network 2). Each line in the event list format contains an event, which consists of the event's start time, end time, first node, and second node, respectively.3. Choose a node similarity file. DynaWAVE accepts node similarities in the node pair list format (see an example node similarity file corresponding to the above two networks). The node similarity file consists of three columns, the first column being nodes from network 1, the second column being nodes from network 2, and the third column being the node similarities. Each node similarity must be between 0 and 1, inclusive. The node similarity file does not necessarily need to contain similarities between all pairs of nodes across the networks.4. Choose the output alignment file name. DynaWAVE will output the alignment file to this file.5. Given network 1 file, ev1.txt, network 2 file, ev2.txt, node similarity file, evsim.txt, and the output alignment file output_alignment.txt, the following command is run in order to align the two networks:julia dynawave.jl ev1.txt ev2.txt evsim.txt output_alignment.txtThe above command will produce the file output_alignment.txt, which contains the alignment. It will also print various statistics related to the alignment."
},

{
    "location": "index.html#secm4-1",
    "page": "Tutorial",
    "title": "IV. Computing node similarities and then running DynaWAVE",
    "category": "section",
    "text": "To follow the instructions below, start julia in the examples/ directory from the base directory of this software, and run the following command to load the DynaWAVE library/module:using DynaWAVE, NetalignUtils, NetalignMeasuresIf you have pre-computed node similarities on which you wish to run DynaWAVE, then you can use the interface at Section III. If you wish to compute the topological node similarities as in the DynaWAVE paper, then perform the instructions in the following sub-section. If you wish to compute node similarities from BLAST E-values to align dynamic protein interaction networks, then go to the sub-section after that."
},

{
    "location": "index.html#sec2-1",
    "page": "Tutorial",
    "title": "Computing topological node similarities",
    "category": "section",
    "text": "Dynamic GDVs [2], i.e., dynamic graphlet degree vectors, are node descriptors that take both the local structural topology and local temporal topology of nodes in a dynamic network into account. The dynamic GDVs of a network can be calculated using the dynamic graphlet counting code available here. In the following example, we will align two networks from the examples/ directory, a dynamic yeast network, and the same network with 10% of its events randomized, while using dynamic GDVs as node similarities.First, we read the networks and the dynamic GDVs into memory.net1 = readeventlist(\"yeastlc_original_tw_1_1.dy\")\ndgdv1 = readgdv(\"yeastlc_original_dgdv_6_4_1.txt\", net1.nodes)\nnet2 = readeventlist(\"yeastlc_rnd_0.10_1_tw_1_1.dy\")\ndgdv2 = readgdv(\"yeastlc_rnd_0.10_1_dgdv_6_4_1.txt\", net2.nodes)Second, we calculate node similarities between node pairs in the two networks using dynamic GDVs, making use of the PCA-based technique described in [2].S = NodeSimMeasure(:pcagdvs,dgdv1,dgdv2).S(Optionally, you can run writelistmat(\"nodesims.txt\", S, net1.nodes, net2.nodes) to save the node similarities to a file, and you can run S = readlistmat(\"nodesims.txt\", net1.nodes, net2.nodes), to read the node similarities from the file.)Finally, we align the two networks using DynaWAVE.f = dynawave(net1.G, net2.G, S)f contains the alignment between the two networks. We can construct the node pairs and calculate node correctness as followsnodepairs = hcat(net1.nodes, net2.nodes[f])\n\nnc = mean(net1.nodes .== net2.nodes[f])We can write the alignment to file as follows.writedlm(\"yeastlc_aln.txt\", nodepairs)\n[2] Exploring the structure and function of temporal networks with dynamic graphlets, Y. Hulovatyy,  H. Chen, and  T. Milenković, Bioinformatics, Volume 31, Issue 12, 15 June 2015, Pages i171–i180, (https://doi.org/10.1093/bioinformatics/btv227)."
},

{
    "location": "index.html#sec3-1",
    "page": "Tutorial",
    "title": "Computing external node similarities from BLAST E-values",
    "category": "section",
    "text": "Instead of dynamic GDVs, we can alternatively use BLAST E-values to align protein interaction networks as follows.First, we read in the E-values to memory. We will be using the same yeast networks as in the above sub-section.E = readevalues(\"yeastlc_yeastlc_evalues.txt\", net1.nodes, net2.nodes)Then, we convert the E-values to node similarities. This is done by converting each E-value to -log(E-value), and then dividing the resulting values with the maximum.S = NodeSimMeasure(:evalues, E).SFinally, we align the two networks.f = dynawave(net1.G, net2.G, S)"
},

{
    "location": "index.html#secm5-1",
    "page": "Tutorial",
    "title": "V. Mimicking experiments from DynaWAVE paper",
    "category": "section",
    "text": "If you wish to create synthetic dynamic networks from network models as in the DynaWAVE paper, then go the following sub-section. If you wish to add noise to a dynamic network as in the DynaWAVE paper, then go to the sub-section after that."
},

{
    "location": "index.html#sec4-1",
    "page": "Tutorial",
    "title": "Creating synthetic dynamic networks",
    "category": "section",
    "text": "In this sub-section, we create random network instances from three network models as in the DynaWAVE paper. Here, we create a random instance of a 1000-node dynamic network using the GEO-GD network model, with parameter p = 03, and linear node arrival. We let the timespan range from 0 to 30 seconds, initializing the network with a 5-node clique.G = rand(GEOGD(0.3, 1, :linear), 1000, 30, 5)Here, we create a random instance of a 1000-node dynamic network using the SF-GD network model, with parameters p = 03 q = 07, and exponential node arrival. We let the timespan range from 0 to 30 seconds, initializing the network with a 5-node clique. G = rand(SFGD(0.3, 0.7, :exp), 1000, 30, 5)Here, we create a random instance of a 1000-node dynamic network using the SNE network model, with parameters lambda = 0032 alpha = 08 beta = 0002, and quadratic node arrival. We let the timespan range from 0 to 30 seconds, initializing the network with a 5-node clique. G = rand(SocialNE(0.032, 0.8, 0.002, :quad), 1000, 30, 5)"
},

{
    "location": "index.html#sec5-1",
    "page": "Tutorial",
    "title": "Adding noise to a dynamic network",
    "category": "section",
    "text": "In this sub-section, given a dynamic network, we add noise to the network. There are two methods we use in the DynaWAVE paper to add noise to a dynamic network: a strict version (strict_events_shuffle) that only changes the event times in the network (from page 15/30 of [3]), and a non-strict version (links_shuffle) that change both event times and links between nodes in the network (from page 16/39 of [3]).Here, we add 30% noise to network G using the strict version.net1 = readeventlist(\"yeastlc_original_tw_1_1.dy\")\nG = net1.G\n\nG30 = strict_events_shuffle(G, 0.30)Here, we add 30% noise to network G using the non-strict version.net1 = readeventlist(\"yeastlc_original_tw_1_1.dy\")\nG = net1.G\n\nG30 = links_shuffle(G, 0.30)[3] Modern temporal network theory: a colloquium, Petter Holme, European Physical Journal B (2015), (https://doi.org/10.1140/epjb/e2015-60657-4)."
},

{
    "location": "static.html#",
    "page": "Static NA",
    "title": "Static NA",
    "category": "page",
    "text": ""
},

{
    "location": "static.html#Static-NA-1",
    "page": "Static NA",
    "title": "Static NA",
    "category": "section",
    "text": ""
},

{
    "location": "static.html#sec6-1",
    "page": "Static NA",
    "title": "Running WAVE with pre-computed node similarities",
    "category": "section",
    "text": "This software tool can also be used to align two static networks using the WAVE method [3]. First, we read in the static networks ex1.gw and ex2.gw (from the examples/ directory). The two networks are very similar, with corresponding nodes of the same names. The networks are in the LEDA format (hence, the .gw extension) and so we use the readgw function to load it.net1 = readgw(\"ex1.gw\")\nnet2 = readgw(\"ex2.gw\")shell> cat ex1.gw\nLEDA.GRAPH\nvoid\nvoid\n-2\n5\n|{A}|\n|{B}|\n|{C}|\n|{D}|\n|{E}|\n5\n1 2 0 |{}|\n2 3 0 |{}|\n3 4 0 |{}|\n4 1 0 |{}|\n4 5 0 |{}|\n\nshell> cat ex2.gw\nLEDA.GRAPH\nvoid\nvoid\n-2\n6\n|{A}|\n|{B}|\n|{C}|\n|{D}|\n|{E}|\n|{F}|\n6\n1 2 0 |{}|\n2 3 0 |{}|\n3 4 0 |{}|\n4 1 0 |{}|\n4 5 0 |{}|\n5 6 0 |{}|Then, we read in the node similarities between the two networks from evgdvsim.txt. The node similarities are also stored in a different format and we use the readdlm function to load it.R = readdlm(\"exgdvsim.txt\", header=true)[1]shell> cat exgdvsim.txt\n5 6\n0.972933 0.93658 0.972933 0.925753 0.935525 0.903174\n0.939805 0.988206 0.939805 0.91356 0.933523 0.937895\n0.972933 0.93658 0.972933 0.925753 0.935525 0.903174\n0.917616 0.925602 0.917616 0.973495 0.910664 0.896031\n0.922085 0.933688 0.922085 0.879783 0.951919 0.954437Finally, we align the two networks using the node similarities.f = wave(net1.G, net2.G, R)Same as before, f contains the alignment between the two networks and we construct the aligned node pairs as follows.nodepairs = hcat(net1.nodes, net2.nodes[f])We can write the alignment to file as follows.writedlm(\"exalnfile.txt\", nodepairs)"
},

{
    "location": "static.html#Computing-node-similarities-and-then-running-WAVE-1",
    "page": "Static NA",
    "title": "Computing node similarities and then running WAVE",
    "category": "section",
    "text": ""
},

{
    "location": "static.html#sec7-1",
    "page": "Static NA",
    "title": "Computing topological node similarities",
    "category": "section",
    "text": "Static GDVs, or GDVs, are node descriptors that take the local topology of nodes in a static network into account. The GDVs of a network can be calculated using GraphCrunch. In the following example, we will align, from the examples/ directory, a yeast network to itself.First, we read the networks and the GDVs into memory.net1 = readgw(\"0Krogan_2007_high.gw\")\ngdv1 = readgdv(\"0Krogan_2007_high.ncount.ndump2\", net1.nodes)\nnet2 = readgw(\"0Krogan_2007_high.gw\")\ngdv2 = readgdv(\"0Krogan_2007_high.ncount.ndump2\", net2.nodes)Second, we calculate node similarities between node pairs in the two networks using GDV similarity.S = NodeSimMeasure(:gdvs, gdv1, gdv2).SFinally, we align the two networks using WAVE.f = wave(net1.G, net2.G, S)We can construct the node pairs and calculate node correctness as followsnodepairs = hcat(net1.nodes, net2.nodes[f])\n\nnc = mean(net1.nodes .== net2.nodes[f])"
},

{
    "location": "static.html#sec8-1",
    "page": "Static NA",
    "title": "Computing external node similarities from BLAST E-values",
    "category": "section",
    "text": "Of course, we can similarly use E-values to align two networks as well. We will be using the same yeast networks as in the above sub-section. We convert the E-values to node similarities and then align.E = readevalues(\"yeastlc_yeastlc_evalues.txt\", net1.nodes, net2.nodes)\nS = NodeSimMeasure(:evalues, E).S\nf = wave(net1.G, net2.G, S)"
},

{
    "location": "static.html#sec9-1",
    "page": "Static NA",
    "title": "Comparison of static and dynamic NA",
    "category": "section",
    "text": "In the DynaWAVE paper, we compared DynaWAVE to WAVE (and MAGNA++ and DynaMAGNA++). The following is an example demonstrating how to make this comparison using WAVE and DynaWAVE.First, we load the two dynamic networks, their corresponding dynamic GDVs, and calculate node similarities.t1 = readeventlist(\"yeastlc_original_tw_1_1.dy\")\ndgdv1 = readgdv(\"yeastlc_original_dgdv_6_4_1.txt\", t1.nodes)\n\nt2 = readeventlist(\"yeastlc_rnd_0.10_1_tw_1_1.dy\")\ndgdv2 = readgdv(\"yeastlc_rnd_0.10_1_dgdv_6_4_1.txt\", t2.nodes)\n\nSdgdv = NodeSimMeasure(:pcagdvs,dgdv1,dgdv2).SThen, we load the corresponding flattened static networks, their corresponding GDVs, and calculate node similarities calculated using the PCA-based technique described in [5].s1 = readgw(\"yeastlc_original.gw\")\ngdv1 = readgdv(\"yeastlc_original.gdv.ndump2\", s1.nodes)\n\ns2 = readgw(\"yeastlc_original.gw\")\ngdv2 = readgdv(\"yeastlc_original.gdv.ndump2\", s2.nodes)\n\nSgdv = NodeSimMeasure(:pcagdvs,gdv1,gdv2).SWe can flatten a dynamic network to a static network as follows. There is an edge between two nodes in the static network if there is atleast one interaction (dynamic edge) between two nodes in the dynamic network.net1 = readeventlist(\"yeastlc_original_tw_1_1.dy\")\nGdynamic = net1.G\n\nGstatic = flatten(Gdynamic)Notice that the flattened version of both t1 and t1 result in the same network as s1 (modulo node permutations) due to our randomization model. That is:flatten(t1[sortperm(t1.nodes)]) == flatten(t2[sortperm(t2.nodes)]) == s1[sortperm(s1.nodes)]\n# trueThen, we align the networks using both WAVE and DynaWAVE. Here, we set the beta parameter to 0.5 as in the paper.fdynamic = dynawave(t1.G, t2.G, Sdgdv, 0.5)\n\nfstatic = wave(s1.G, s2.G, Sgdv, 0.5)We can calculate node correctness as follows.nc_dynamic = nodecorrectness(fdynamic, t1.nodes, t2.nodes)\n\nnc_static = nodecorrectness(fstatic, s1.nodes, s2.nodes)Finally, above, we aligned the yeast dynamic network to the same network with 10% randomization. We can also align the yeast dynamic network to itself.S = NodeSimMeasure(:pcagdvs,gdv1,gdv1).S\nf = dynawave(t1.G, t1.G, S, 0.5)\nnodecorrectness(f, t1.nodes, t1.nodes)"
},

{
    "location": "funs.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "funs.html#DynaWAVE.dynawave",
    "page": "Functions",
    "title": "DynaWAVE.dynawave",
    "category": "Function",
    "text": "dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,\n     S::AbstractMatrix, [ beta=0.5,seeds=[] ]);\n     skipalign=false,details=false) -> f [, M]\n\nGiven two dynamic networks and node similarities between them, align the two networks by running the DynaWAVE algorithm.\n\nArguments\n\nG1,G2 : input networks in sparse matrix format\nS : node similarities beween the two networks\nbeta : weight between edge and node conservation, beta=1.0 weighs   edge conservation highly while beta=0.0 weighs node conservation highly\nseeds : seed aligned node pairs. For example, if we know that the   3rd node in the first network is aligned to the 7th node in the   second network, and the 5th node in the first network is aligned to the   9th node in the second network, then we will set seeds = [(3,7), (5,9)]\n\nKeyword arguments\n\nskipalign : Don't align; just return the alignment details in WaveModel.\ndetails : Returns M, the WaveModel.\n\nOutput\n\nf : Alignment, i.e. node mapping from G1 to G2. f[i] describes   node pairnodes1[i], nodes2[f[i]], wherenodes1andnodes2are   the vectors of node names corresponding toG1andG2` respectively.\n\n\n\n"
},

{
    "location": "funs.html#DynaWAVE.wave",
    "page": "Functions",
    "title": "DynaWAVE.wave",
    "category": "Function",
    "text": "wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,\n     S::AbstractMatrix, [ beta=0.5,seeds=[] ]);\n     skipalign=false,details=false) -> f [, M]\n\nGiven two static networks and node similarities between them, align the two networks by running the WAVE algorithm (Yihan Sun, Joseph Crawford, Jie Tang, and Tijana Milenkovic, Simultaneous Optimization of Both Node and Edge Conservation in Network Alignment via WAVE, in Proceedings of the Workshop on Algorithms in Bioinformatics (WABI), Atlanta, GA, USA, September 10-12, 2015, pages 16-39).\n\nArguments\n\nG1,G2 : Input networks in sparse matrix format.\nS : Node similarities beween the two networks.\nbeta : Weighs between edge and node conservation, beta=1.0 weighs   edge conservation highly while beta=0.0 weighs node conservation highly\nseeds : Seed aligned node pairs. For example, if we know that the   3rd node in the first network is aligned to the 7th node in the   second network, and the 5th node in the first network is aligned to the   9th node in the second network, then we will set seeds = [(3,7), (5,9)].\n\nKeyword arguments\n\nskipalign : Don't align; just return the alignment details in WaveModel.\ndetails : Returns M, the WaveModel.\n\nOutput\n\nf : Alignment, i.e. node mapping from G1 to G2. f[i] describes   node pairnodes1[i], nodes2[f[i]], wherenodes1andnodes2are   the vectors of node names corresponding toG1andG2` respectively.\n\n\n\n"
},

{
    "location": "funs.html#DynaWAVE.shufflealign",
    "page": "Functions",
    "title": "DynaWAVE.shufflealign",
    "category": "Function",
    "text": "shufflealign(method, G1::SparseMatrixCSC,G2::SparseMatrixCSC,\n                  S::AbstractMatrix, beta::Float64,\n                  seeds=[], details=false) -> f [, M]\n\nGiven two networks and node similarities between them, before aligning using either dynawave or wave, this shuffles each input network independently before passing the networks to dynawave or wave. This is to remove node order biases when the two networks being aligned have similar node sets with similar node orderings. This is for use in evaluation, where biases like this tend to show much better results that is actually possible with a particular method. We obviously do not want biases like these during evaluation.\n\nFor example, when aligning a network to itself evaluate a network alignment method, the node order of the two networks will be the same and the node set will be the same. This results in the network alignment method \"knowing\" the true node mapping, and often producing alignments of much higher quality that is actually possible with the method were the true node mapping not known. This method randomizes (shuffles) the node order before handing it to dynawave or wave and then returns the deshuffled alignment.\n\nArguments\n\nmethod : If method = dynawave, run the DynaWAVE algorithm after   shuffling the networks. If method = wave, run the WAVE algorithm   after shuffling the networks.\nSee dynawave and wave for the other arguments.\n\nOutput\n\nf : The resulting alignment has the correct node order with respect to the input G1 and G2 networks; i.e., it is \"deshuffled\".\n\n\n\n"
},

{
    "location": "funs.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "CurrentModule = DynaWAVEdynawave\nwave\nshufflealign"
},

{
    "location": "internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals.html#DynaWAVE.WaveModel",
    "page": "Internals",
    "title": "DynaWAVE.WaveModel",
    "category": "Type",
    "text": "WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,\n          meas::NetalignMeasure, seeds = [])\n\nInternal function that creates a problem structure given two networks (static or dynamic) to align by and NetalignMeasure to optimize over. See align! on how perform the alignment.\n\nArguments\n\nG1,G2 : input networks in sparse matrix format\nmeas : NetalignMeasure object (see NetalignMeasures.jl). WaveModel can optimize the DWEC measure (DWECMeasure), a dynamic network alignment measure and the WEC measure (WECMeasure), a static network alignment method.\n\n\n\n"
},

{
    "location": "internals.html#DynaWAVE.align!-Tuple{DynaWAVE.WaveModel}",
    "page": "Internals",
    "title": "DynaWAVE.align!",
    "category": "Method",
    "text": "align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-length(M.seeds)),\n            verbose=false) -> f\n\nInternal function that performs alignment with respect to the parameters in WaveModel.\n\nArguments\n\nmaxiter : Stop aligning after maxiter iterations. Each iteration   creates one aligned node pair.\nverbose : Prints debugging outputs\n\n\n\n"
},

{
    "location": "internals.html#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": "CurrentModule = DynaWAVEWaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,\n                   meas::NetalignMeasure, seeds = zeros(Int,(0,0)))\nalign!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-size(M.seeds,1)),\n                verbose=false,details=false)"
},

]}
