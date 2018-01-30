# DynaWAVE

[![Build
 Status](https://travis-ci.org/vvjn/DynaWAVE.jl.svg?branch=master)](https://travis-ci.org/vvjn/DynaWAVE.jl)
 [![codecov.io](http://codecov.io/github/vvjn/DynaWAVE.jl/coverage.svg?branch=master)](http://codecov.io/github/vvjn/DynaWAVE.jl?branch=master)  [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://vvjn.github.io/DynaWAVE.jl/latest)


Alignment of dynamic (temporal) networks with the DynaWAVE
method. Network alignment aims to find a mapping between networks that
conserves similar regions between the networks. Dynamic (or temporal)
networks are networks whose interactions are time dependent. DynaWAVE
aligns two dynamic networks while taking temporal information into
account. An implementation of the WAVE method, which aligns two static
networks is also available.

Home page: https://www.nd.edu/~cone/DynaWAVE/ Reference: Vipin Vijayan
and Tijana Milenković, Aligning dynamic networks with DynaWAVE,
Bioinformatics,
[btx841](https://doi.org/10.1093/bioinformatics/btx841), 2017.

# Installation

DynaWAVE can be installed by starting julia and executing the
following command. This will install DynaWAVE and packages that it
depends on.

```julia
Pkg.add("DynaWAVE")
```

# Documentation

Available [here](https://vvjn.github.io/DynaWAVE.jl/latest).
