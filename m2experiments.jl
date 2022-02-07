R = QQ[x, y, t]
f = 10^5 * t^2 + 10^5*t + 100
g = 1/100 * t^2 + 100 * t + 10
I = ideal(x - f, y - g)
eliminate(I, t)
F = x^2  - 20000000x*y + 100000000000000y^2   -9999800000200x - 1999998000000000y + 10999980000010000
expand(subs(F, [x, y] => [f, g]))

R = QQ[x, y, t]
f = 10^2 * t^2 + 10^(-2) * t + 10^(-7)
g = 10^(4) * t^2 - t + 10^(-5)
I = ideal(x - f, y - g)
eliminate(I, t)



ideal(x^2  - 20000000x*y + 100000000000000y^2 
- 9998800000200x +7999002000000000y - 88990119999990000)


function xx(a)
     maximum([5+2*a;5+a;2])
end
function yy(a)
     maximum([-2+2*a;2+a;1])
end


as = -10:0.01:10


plot!(xx.(as),yy.(as))


function ϕ(a)
     [10. ^(5+2*a).+10. ^(5+a).+10^2; (10).^(-2. +2*a).+10. ^(2. +a).+10^1]
end
function ψ(a)
     r = 10*randn(ComplexF64)
     [10^5*r^(2*a)+10^5*r^a+10^2;10^(-2)*r^(2*a)+10^2*r^a+10^1]
end

vals2 = log10.(abs.((transpose(hcat(ψ.(as)...)))));
scatter!(vals2[:,1],vals2[:,2],markersize = 1.1, color = "blue")




pt = [xx.(as)[1900],yy.(as)[1900]]
affine_lin_fcts =
[(log10(1), [2,0]),
(log10(20000000) -0.3, [1,1]),
(log10(100000000000000), [0,2]),
(log10(9998980000200), [1,0]),
(log10(7999002000000000), [0,1]),
(log10(88990119999990000), [0, 0])]






function sample_from_trop(affine_lin_fcts, n)
    pts = []
    for i in 1:n
        a, b, aa, bb = BigFloat.(randn(Float64, 4))
        for i in 1:length(affine_lin_fcts), j in 1:length(affine_lin_fcts)
            f, g = affine_lin_fcts[[i,j]]
            if i != j
                t = ( g[1] - f[1] + [b, bb]'*(g[2] - f[2]) )/([a, aa]'*(f[2] - g[2]))
                vt = [a*t + b, aa*t + bb]
                if 0 < vt[1] && vt[1] < 20 && -10 < vt[2] && vt[2] < 20
                    value, ind = findmax([f[1] + vt'f[2] for f in affine_lin_fcts])
                    #println("t = $t")
                    #println("v(t) = $vt")

                    #println("optimizer = ", affine_lin_fcts[ind])
                    #pts = push!(pts, vt)
                    if ind in [i, j]
                        pts = push!(pts, vt)
                    end
                end
            end
        end
    end
    x = [pt[1] for  pt in pts]
    y = [pt[2] for  pt in pts]
    x, y
end

v1 = [log10(BigFloat(1)),2,0]
v2 = [log10(BigFloat(20000000)),1,1]
v3 = [log10(BigFloat(100000000000000)),0,2]

n = 5000
x, y = sample_from_trop(affine_lin_fcts, 10000)
scatter(x, y, markersize = 1.1, markercolor = :green, markerstrokewidth = 0.5)
scatter!(vals2[:,1],vals2[:,2],markersize = 1.1, color = "blue")
plot!(xx.(as),yy.(as))







@polyvar x y
f = x^2*y + exp(1)*y^2 + exp(4)*x^2*y^3 + exp(-4)*y^4 + exp(2)*x^4*y^4
A = amoeba(f);
plot(A)
amoeba_spine_plot = plot!(spine(f))


f = x^2  - 20000000*x*y + 100000000000000*y^2  - 9998800000200*x +
 9799002000000000*y - 8899101999990000

