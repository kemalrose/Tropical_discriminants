

R = QQ[x, y, t]
f = (t - 1/25)^3 * (t - 5^3)
g = (t - 1/25)^1 * (t - 5)^2  * (t - 5^3)
I = ideal(x * (t - 5)^3 * (t - 25)^3 - f, y * (t - 25)^4- g)
toString eliminate(I, t)
F = 1490116119384765625000000000000*x^4*y^3-389074785232543945312500000000*x^3*y^3+4106795716054687500000000000*x^2*y^4-6031830555106200000000000*x*y^5+2137378418265692160000*y^6+27119010589599609375000000000*x^3*y^2+15783626283118542480468750000*x^2*y^3+37602347079547355625000000*x*y^4-4274756836531384320000*y^5+2140639939707144375000000000*x^2*y^2-34191816324780610291984375*x*y^3+2137378418265692160000*y^4+47403977291875310625000000*x^2*y-2960714274258647127250*x*y^2+256993127401722120000000*x^2-16127138372745043831*x*y


expand(subs(F, [x, y] => [f, g]))






R = QQ[a, b, c, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9]
x0 = a^3 - t0
x1 = a^2 * b - t1
x2 = a^2 * c - t2
x3 = a * b^2 - t3
x4 = a * b * c - t4
x5 = a * c^2 - t5
x6 = b^2 * c - t6
x7 = b^3 - t7
x8 = b * c^2 - t8
x9 = c^3 - t9
I = ideal(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9)
eliminate({a, b, c},I)



R = QQ[a, b, c, t0, t2, t3, t6, t8, t9]
x0 = a^3 - t0
x2 = a^2 * c - t2
x3 = a * b^2 - t3
x6 = b^2 * c - t6
x8 = b * c^2 - t8
x9 = c^3 - t9
I = ideal(x0, x2, x3, x6, x8, x9)
eliminate({a, b, c},I)



I = ideal(x - f, y - g)
eliminate(I, t)
F = x^2  - 20000000x*y + 100000000000000y^2   -9999800000200x - 1999998000000000y + 10999980000010000
expand(subs(F, [x, y] => [f, g]))







R = QQ[u0, u1, x1, x2, x3, x4]
I = ideal(
x1 - (u1 - u0 * 5^(-2)),
x2 - (u1 - u0 * 5^(1)),
x3 - (u1 - u0 * 5^(3)),
x4 - u0
)
J = ideal(
x1 - (u1 - u0 * 5^(-2)),
x2 - (u1 - u0 * 5^(1)),
x3 - (u1 - u0 * 5^(3)),
x4 - u0,
u0 - 1
)
toString gfanTropicalBasis I
  {(1/25)*u0-u1+x1, 5*u0-u1+x2, 125*u0-u1+x3, -u0+x4,
      u0-(1/120)*x2+(1/120)*x3, u0-(25/124)*x1+(25/124)*x2,
      u0-(25/3124)*x1+(25/3124)*x3}



R = QQ[x, y, t]
f = (t - 5^(-2))^2
g = (t - 5^(-2))^1 * (t - 5^(3))^(1) 
I = ideal(x * (t - 5^(3))^(1) * (t - 5^(1)) - f, y * (t - 5^(1))^(2) - g )
eliminate(I, t)
152587890625*x^4*y^2-610351562500*x^3*y^2-35728465575000000000*x^3*y
      +915527343750*x^2*y^2-7444919166800000000*x^2*y-610351562500*x*y^2-
      300066291800000000*x*y+152587890625*y^2-893802369030533349376*x


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


function ??(a)
     [10. ^(5+2*a).+10. ^(5+a).+10^2; (10).^(-2. +2*a).+10. ^(2. +a).+10^1]
end
function ??(a)
     r = 10*randn(ComplexF64)
     [10^5*r^(2*a)+10^5*r^a+10^2;10^(-2)*r^(2*a)+10^2*r^a+10^1]
end

vals2 = log10.(abs.((transpose(hcat(??.(as)...)))));
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






 R = QQ[p]
 F = R[x]

