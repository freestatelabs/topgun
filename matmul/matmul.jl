using BenchmarkTools

function mul!(c, a, b, m, n, p) 
    # C, A, B are all column-major 
    #
    # To access an element of a column-major array stored in a 1d vector, 
    # A[i,j] = a[(j-1)*row_length+i]

    for i in 1:m 
        for j in 1:p
            for k in 1:n
                c[(j-1)*p+i] += a[(k-1)*n+i]*b[(j-1)*p+k]
            end
        end
    end

end


function mul2!(C, A, B, m, n, p) 

    for i in 1:m 
        for j in 1:p 
            for k in 1:n 
                C[i,j] += A[i,k]*B[k,j]
            end
        end
    end
end


m = 100
n = 100
p = 100
a = rand(m*n)
b = rand(n*p)
c = zeros(m*p)

mul!(c,a,b,m,n,p)

# x = @benchmark mul!($c,$a,$b,$m,$n,$p)
# julia> @benchmark mul!($c,$a,$b,$m,$n,$p)
# BenchmarkTools.Trial: 6626 samples with 1 evaluation per sample.
#  Range (min … max):  710.666 μs …   5.665 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     725.125 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   752.811 μs ± 176.278 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

#   █▅▅▅▅▄▃▃▂▂▂▂▁▁▁▁▁▁ ▁                                          ▁
#   █████████████████████████▇▇▇▆▇▆▆▆▇▆▆▇▅▆▅▆▅▅▅▅▅▄▅▄▅▅▃▅▅▁▄▃▃▄▃▅ █
#   711 μs        Histogram: log(frequency) by time       1.08 ms <

#  Memory estimate: 0 bytes, allocs estimate: 0.


A = zeros(m,n) 
B = zeros(n,p) 
C = zeros(m,p)

for i in 1:m*n 
    A[i] = a[i] 
end

for i in 1:n*p 
    B[i] = b[i] 
end 

function julia_mul!(C, A, B) 
    C .= A*B 
end 

# y = @benchmark julia_mul!(C, A, B)
# julia> @benchmark julia_mul!(C, A, B)
# BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
#  Range (min … max):  51.083 μs …  6.034 ms  ┊ GC (min … max): 0.00% … 98.37%
#  Time  (median):     60.917 μs              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   66.098 μs ± 99.447 μs  ┊ GC (mean ± σ):  5.99% ±  5.97%

#          ▅█▅▁                                                  
#   ▃▄▅▅▄▄▇████▆▅▄▃▃▃▂▂▂▂▂▂▂▂▂▂▁▂▂▁▂▂▂▁▂▂▁▂▂▂▂▂▁▁▁▂▂▂▂▁▁▂▁▁▁▁▁▂ ▃
#   51.1 μs         Histogram: frequency by time         119 μs <

#  Memory estimate: 78.20 KiB, allocs estimate: 3.


julia_mul!(C, A, B)

print("Error: ")
function checkerror(C, c, m, p)
    maxerr = 0.0 
    for i in 1:m*p 
        err = abs(C[i] - c[i])
        if err > maxerr
            maxerr = err 
        end 
    end 
    return maxerr
end
print(checkerror(C, c, m, p))
println("")
