using Plots, BenchmarkTools, Printf, Statistics, LaTeXStrings

function ellipK(k2::Real; itmax=100, errmax=1e-12)


    it = 0
    err = 2*errmax

    a0 = 1.0; g0 = sqrt(1-k2)

    while it < itmax && err > errmax

        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0)
        a0 = a1; g0 = g1

        err = abs(a0-g0)
        it += 1

    end

    return pi/(2*a0), it
end

function ellipK2(k2::Real)

    a0 = 1.0; g0 = sqrt(1-k2)
    if k2 < 0.94
        if k2 < 0.33 
            if k2 < 0.01
                # 2 iterations 
                a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
                a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
            else
                # 3 iterations 
                a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
                a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
                a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
            end
        else
            # 4 iterations 
            a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
            a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
            a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
            a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
        end

    else 
        # 6 iterations 
        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
        a1 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
    end

    return pi/(2*a0), 0
end


function ellipK3(k2::Real)

    a0 = 1.0; g0 = sqrt(1-k2)
    a0 = 0.5*(a0+g0); g0 = sqrt(a0*g0); a0 = a1; g0 = g1
    a0 = 0.5*(a0+g0); g0 = sqrt(a0*g0); a0 = a1; g0 = g1
    a0 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
    a0 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
    a0 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1
    a0 = 0.5*(a0+g0); g1= sqrt(a0*g0); a0 = a1; g0 = g1

    return pi/(2*a0)
end


function ellipE(k2::Real; itmax=100, errmax=1e-8)

    err = 2*errmax
    n = 0
    an = 1; gn = sqrt(1-k2); cn = abs(an^2 - gn^2)
    esumn = cn*2.0^(n-1)

    while n < itmax && err > errmax
        n +=1 

        an1 = (an + gn)/2.0
        gn1 = sqrt(an*gn)
        cn1 = abs(an1^2 - gn1^2)
        esumn1 = esumn + cn1*(2^(n-1))

        err = abs(esumn1 - esumn)

        an = an1 
        gn = gn1 
        esumn = esumn1 

    end

    return (1-esumn)*pi/(2*an), n
end


N = 10000
K = zeros(N); K2 = zeros(N); E = zeros(N); E2 = zeros(N)
Kit = zeros(N); Kit2 = zeros(N); Eit = zeros(N) 
k2 = LinRange(0, 1.0-1e-12, N)

for i = 1:N 
    K[i], Kit[i] = ellipK(k2[i])
    K2[i], _ = ellipK2(k2[i])
    E[i], Eit[i] = ellipE(k2[i])
    # @printf "k2 = %f, Kit = %i, Eit = %i\n" k2[i] Kit[i] Eit[i]
end

p = plot(k2, Kit, label=L"K(k^2)")
plot!(p, k2, Eit, label=L"E(k^2)")
xlabel!(L"k^2")
ylabel!(L"Iterations")
title!("Convergence Criteria for Complete Elliptic Integrals")
display(p)
savefig("convergence.svg")
# @benchmark ellipK3.($k2)
# @benchmark ellipE.($k2)

@printf "Max error in K and K2: %f\n" maximum(K .- K2)
@printf "Avg error in K and K2: %f\n" mean(K .- K2)

# Breakpoints for K 
# 2/3: 0.01 
# 3/4: 0.33 
# 4/5: 0.94


