# Internals

```@meta
CurrentModule = DynaWAVE
```

``` @docs
WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                   meas::NetalignMeasure, seeds = zeros(Int,(0,0)))
align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-size(M.seeds,1)),
                verbose=false,details=false)
```
