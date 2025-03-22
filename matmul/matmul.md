# Matrix Multiplication

Head to head compiler showdown. Test matrix multiplication in different languages

## Languages
- C 
- C++
- Fortran
- Rust
- Julia

## Algorithm

``` 
C = A*B

Given:
A: m x n matrix 
B: n x p matrix 
C: m x p matrix 

for i in 1:m 
    for j in 1:p 
        for k in 1:n 
            C[i,j] += A[i,k]*B[k,j] 
```
## Permutations
1. Naive implementation 
2. Implementation using language features
