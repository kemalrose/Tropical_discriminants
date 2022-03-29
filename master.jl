using Arpack, Oscar, HomotopyContinuation, LinearAlgebra
import Pkg; Pkg.add("MultiFloats")
using GenericSVD
using Plots, PolynomialAmoebas

MultiFloats.use_bigfloat_transcendentals()
#Pkg.add("https://github.com/JuliaLinearAlgebra/GenericSVD.jl.git")




include("get_vtx.jl")
include("get_newton.jl")
include("tropicalize_curve.jl")
include("plot_trop_curve.jl")


p = 5
θ = [(5//1)^(-2), 5, 5^2, 5^3]
m = [3, -3, -3, 1]
n = [1, 2, -4, 1]

verts, max_pols  = trop_curve(θ, m, n, p)

plt = plot_curve(verts, max_pols)
plot(plt, legend = false)





R = QQ[x, y, t]
ϕ1 = (t - 1/25)^3 * (t - 5^3)
ϕ2 = (t - 1/25)^1 * (t - 5)^2  * (t - 5^3)
I = ideal(x * (t - 5)^3 * (t - 25)^3 - ϕ1, y * (t - 25)^4- ϕ2)
toString eliminate(I, t)
f = 1490116119384765625000000000000*x^4*y^3-389074785232543945312500000000*x^3*y^3+4106795716054687500000000000*x^2*y^4-6031830555106200000000000*x*y^5+2137378418265692160000*y^6+27119010589599609375000000000*x^3*y^2+15783626283118542480468750000*x^2*y^3+37602347079547355625000000*x*y^4-4274756836531384320000*y^5+2140639939707144375000000000*x^2*y^2-34191816324780610291984375*x*y^3+2137378418265692160000*y^4+47403977291875310625000000*x^2*y-2960714274258647127250*x*y^2+256993127401722120000000*x^2-16127138372745043831*x*y




@polyvar x y t
mons = DynamicPolynomials.monomials((1+x)^10*(1+y)^8)
f = dot(exp.(10*rand(length(mons))),mons)
f = 1490116119384765625000000000000*x^4*y^3-389074785232543945312500000000*x^3*y^3+4106795716054687500000000000*x^2*y^4-6031830555106200000000000*x*y^5+2137378418265692160000*y^6+27119010589599609375000000000*x^3*y^2+15783626283118542480468750000*x^2*y^3+37602347079547355625000000*x*y^4-4274756836531384320000*y^5+2140639939707144375000000000*x^2*y^2-34191816324780610291984375*x*y^3+2137378418265692160000*y^4+47403977291875310625000000*x^2*y-2960714274258647127250*x*y^2+256993127401722120000000*x^2-16127138372745043831*x*y


mons = DynamicPolynomials.monomials(f)
exps = exponents.(DynamicPolynomials.terms(f))
coeffs = DynamicPolynomials.coefficients.(f)
vtcs = hcat(exps...)'
#weights = rationalize.(log10.(abs.(coeffs)))
QQ5 = FlintPadicField(5,20)

weights = valuation.(QQ5.(coeffs)).*(-1)
p, maxcells, vertsoftrop, upperfacets = plot2Dsubdivision(vtcs,weights)
plot(p)

p = plotTropicalCurve(vertsoftrop,vtcs,upperfacets,maxcells)

plotUpperHull(vtcs, weights)