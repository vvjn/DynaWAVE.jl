# Functions

```@meta
CurrentModule = DynaWAVE
```

``` @docs
dynawave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                  S::AbstractMatrix, beta=0.5,
                  seeds=zeros(Int,(0,0));skipalign=false,details=false)
wave(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
              S::AbstractMatrix, beta=0.5,
              seeds=zeros(Int,(0,0));skipalign=false,details=false)
shufflealign(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                      S::AbstractMatrix, beta::Float64,
                      seeds=Vector{Tuple{Int,Int}}();
                      details=false,method=dynawave)
WaveModel(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                   meas::NetalignMeasure, seeds = zeros(Int,(0,0)))
align!(M::WaveModel; maxiter::Integer=(size(M.G1,1)-size(M.seeds,1)),
                verbose=false,details=false)
```
