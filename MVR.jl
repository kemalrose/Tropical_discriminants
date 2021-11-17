using HomotopyContinuation

X = 2
XX = 2
M = 2
MM = 2
N = 4.0

@var A[1:M,1:X,1:MM,1:XX] μ[1:X] ν[1:XX] ϕ[1:X] ψ[1:XX]

eqns = fill(zero(sum(A) + sum(μ) + sum(ϕ) + sum(ψ)),X + XX)

for x = 1:X
    for m = 1:M
        for mm = 1:MM
            for xx = 1:XX
                eqns[x] += mm*ϕ[x]^m*A[m,x,mm,xx]*ψ[xx]^mm
            end
        end
    end
    eqns[x] += - μ[x]
end

for xx = 1:XX
    for m = 1:M
        for mm = 1:MM
            for x = 1:X
                eqns[X + xx] += m*ϕ[x]^m*A[m,x,mm,xx]*ψ[xx]^mm
            end
        end
    end
    eqns[X + xx] += - ν[xx]
end

parvec = [A[:];μ;ν]

Sys = System(eqns,parameters = parvec)

R = solve(Sys;target_parameters = rand(length(parvec)))

function get_real_solutions(R)
    realsols = real_solutions(R)
    posinds = findall(sol -> prod(sol .> 0) == 1, realsols)
    return realsols[posinds]
end
