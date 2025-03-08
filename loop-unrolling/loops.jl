

using BenchmarkTools

function f1(x) 

    y = 0 
    for _x in x 
        y += _x
    end 
    return y
end

function f2(x; skip=2) 
    y = 0
    Nmax = length(x) - (length(x) % (skip-1))
    _x = zeros(skip)
    for i in 1:skip:Nmax

        for j in 0:(skip-1)
            _x[j+1] = x[i+j] 
        end 


        y += sum(_x)
    end

    if (Nmax % skip) > 0 
        for i in (Nmax+1):(Nmax + (Nmax % skip) )
            y += x[i] 
        end 
    end 
    return y
end


x = rand(1_000_002)

# t = @benchmark f1($x) 
# display(t)
# t = @benchmark f2($x) 
# display(t)

f1(x)
f2(x)