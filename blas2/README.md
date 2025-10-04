# BLAS2

Tests different variants of matrix-vector multiply. This has been a bottleneck on `shepherd`.

## References
https://discourse.julialang.org/t/why-mul-is-so-fast/31436/6
https://github.com/flame/how-to-optimize-gemm/wiki#the-gotoblasblis-approach-to-optimizing-matrix-matrix-multiplication---step-by-step

## Reference Solution
Reference solution is julia 1.11.7: 
```julia
using BenchmarkTools, LinearAlgebra

n = 10000
A = rand(n, n)
x = rand(n)
b = rand(n)

res = @benchmark mul!($b, $A, $x)
display(res)
```
Result:
```bash
BenchmarkTools.Trial: 391 samples with 1 evaluation per sample.
 Range (min … max):  10.722 ms … 24.821 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     12.574 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   12.748 ms ±  1.187 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

              ▂█         ▁▁                                    
  ▄▃▂▃▃▃▃▃▃▄▃▄██▄▆▆▇▇▅▃▅▇██▇▄▂▃▂▄▇▆▂▂▁▂▁▁▃▁▁▂▁▁▁▂▁▁▁▂▂▁▁▁▁▂▁▂ ▃
  10.7 ms         Histogram: frequency by time        16.7 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.
 ```

Oops, forgot to turn off BLAS multithreading...
```julia
BLAS.set_num_threads(1)
```

BenchmarkTools.Trial: 152 samples with 1 evaluation per sample.
 Range (min … max):  32.542 ms … 43.190 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     32.695 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   32.942 ms ±  1.102 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

    ▅█▂▂                                                       
  ▃▅████▅▃▂▄▄▃▄▂▃▂▃▂▃▂▁▁▄▃▃▁▂▂▂▂▃▁▂▁▁▁▁▁▁▂▁▁▁▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▂ ▂
  32.5 ms         Histogram: frequency by time        34.4 ms <

 Memory estimate: 0 bytes, allocs estimate: 0.

 So, **33ms** is the time to beat. 

 ## Baseline Solution
 Nested arrays for the matrix, i.e.:

 ```c
 // A*x = b
int mmultv(Matrix *A, Vector *x, Vector *b) {
    if (A->ncols != x->nrows) { return -1; }
    for (size_t i=0; i<A->nrows; i++) {
        double acc = 0.0;
        for (size_t j=0; j<A->ncols; j++) {
            acc += A->values[i][j]*x->values[j];
            //printf("i, j = %li, %li\n", i, j);
        }
        b->values[i] = acc;
    }
    return 0;
}
```

Speed: 122.9 ms

## Linear Array Solution 
Access A via A[row][col] = A[row * ncols + col]:
```c
int mmultv(Matrix *A, Vector *x, Vector *b) {
    if (A->ncols != x->nrows) { return -1; }
    for (size_t i=0; i<A->nrows; i++) {
        double acc = 0.0;
        for (size_t j=0; j<A->ncols; j++) {
            acc += A->values[i*A->ncols + j]*x->values[j];
        }
        b->values[i] = acc;
    }
    return 0;
}
```

Speed: still 123 ms