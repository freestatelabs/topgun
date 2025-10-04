using BenchmarkTools, LinearAlgebra
BLAS.set_num_threads(4)

n = 10000
A = rand(n, n)
x = rand(n)
b = rand(n)

res = @benchmark mul!($b, $A, $x)
display(res)
