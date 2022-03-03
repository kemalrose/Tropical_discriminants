using Oscar
using DynamicPolynomials
using Plots, Statistics

function orderVtcs(vtcs)
    meanvert = mean(vtcs,dims=1)
    diffvecs = vtcs.-meanvert
    per = sortperm([diffvecs[i,2]/diffvecs[i,1] for i ∈ 1:size(diffvecs,1)])
    return vtcs[per,:]
end

function plot2Dsubdivision(vtcs,weights)
    lifted_points = [[vtcs[i,:];weights[i]] for i = 1:size(vtcs,1)]
    lifted_points = hcat(lifted_points...)'
    P = convex_hull(lifted_points)
    Σ = normal_fan(P)
    R = rays(Σ)
    maxcells = []
    vertsoftrop = []
    upperfacets = []
    for i = 1:length(R)
        if R[i][end] < 0
            w = [dot(R[i],lifted_points[j,:]) for j = 1:size(lifted_points,1)]
            vertsincell = findall(j->w[j]==minimum(w),1:length(w))
            push!(maxcells,vertsincell)
            push!(vertsoftrop, R[i][1:end-1].//R[i][end])
            push!(upperfacets,convex_hull(lifted_points[vertsincell,:]))
        end
    end
    i = 1 
    v = vtcs[maxcells[i],:]
    v = convert.(Int64,hcat(collect(vertices(convex_hull(v)))...)')
    v = orderVtcs(v)
    p = plot([v[:,1]; v[1,1]],[v[:,2];v[1,2]],color = "black",legend = false)
    for i = 2:length(maxcells)
        v = vtcs[maxcells[i],:]
        plot!(p, [v[:,1]; v[1,1]],[v[:,2];v[1,2]],color = "black")
    end
    return p, maxcells, vertsoftrop, upperfacets
end

function plotUpperHull(vtcs,weights; mycolor = "black")
    p, maxcells = plot2Dsubdivision(vtcs,weights)
    i = 1
    v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
    v = convert.(Rational{BigInt},hcat(collect(vertices(convex_hull(v)))...)')
    v = orderVtcs(v)
    q = plot([v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = mycolor,legend = false)
    plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = "blue")
    for i = 1:length(maxcells)
        v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
        v = convert.(Rational{BigInt},hcat(collect(vertices(convex_hull(v)))...)')
        v = orderVtcs(v)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = mycolor,legend = false)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = "blue")
    end
    return q
end

function plotTropicalCurve(vertsoftrop,vtcs,upperfacets,maxcells; margin = 1.0)
    xmin = minimum([v[1] for v ∈ vertsoftrop]) - margin
    xmax = maximum([v[1] for v ∈ vertsoftrop]) + margin
    ymin = minimum([v[2] for v ∈ vertsoftrop]) - margin
    ymax = maximum([v[2] for v ∈ vertsoftrop]) + margin
    m = length(vertsoftrop)
    edgesoftrop = []

    if length(upperfacets) == 1
        p = plot()
    else
        for i = 1:m
            for j = i+1:m
                if dim(intersect(upperfacets[i],upperfacets[j])) == 1
                    push!(edgesoftrop,[i j])
                end
            end
        end
    e = edgesoftrop[1]
    H = hcat([Float64.(vertsoftrop[e[i]]) for i = 1:2]...)
    p = plot(H[1,:],H[2,:], color = "black", legend = false)
    for e ∈ edgesoftrop[2:end]
        H = hcat([Float64.(vertsoftrop[e[i]]) for i = 1:2]...)
        plot!(p,H[1,:],H[2,:], color = "black")
    end
    end
    P = convex_hull(vtcs)
    edgesofP = faces(P,1)
    verticesofP = vertices(P)
    for e ∈ edgesofP
        for i =1:length(maxcells)
            if dim(intersect(e,convex_hull(vtcs[maxcells[i],:]))) == 1
                v = vertices(e)
                w = v[1]-v[2]
                taildirection = [-w[2];w[1]]
                if dot(taildirection,v[1]) != maximum([dot(taildirection,vv) for vv ∈ verticesofP])
                    taildirection = -taildirection
                end
                if taildirection[1] > 0 
                    k1 = (xmax - vertsoftrop[i][1])/taildirection[1]
                elseif taildirection[1]<0
                    k1 = (xmin - vertsoftrop[i][1])/taildirection[1]
                else
                    k1 = Inf
                end
                if taildirection[2] > 0 
                    k2 = (ymax - vertsoftrop[i][2])/taildirection[2]
                elseif taildirection[2] < 0
                    k2 = (ymin - vertsoftrop[i][2])/taildirection[2]
                else
                    k2 = Inf
                end
                q = hcat(vertsoftrop[i],vertsoftrop[i]+minimum([k1;k2])*taildirection)
                plot!(p,q[1,:],q[2,:],color = "black")
            end
        end
    end
    return p
end



function plotUpperHull(vtcs,weights; hullcolor = "black", groundcolor = "blue")
    p, maxcells = plot2Dsubdivision(vtcs,weights)
    i = 1
    v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
    v = convert.(Rational{BigInt},hcat(collect(vertices(convex_hull(v)))...)')
    v = orderVtcs(v)
    q = plot([v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = hullcolor,legend = false)
    plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = groundcolor)
    for i = 1:length(maxcells)
        v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
        v = convert.(Rational{BigInt},hcat(collect(vertices(convex_hull(v)))...)')
        v = orderVtcs(v)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = hullcolor,legend = false)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = groundcolor)
    end
    return q
end



@polyvar x y t


R = QQ[x, y, t]
f = (t - 5^(-2))^2
g = (t - 5^(-2))^1 * (t - 5^(3))^(1) 
I = ideal(x * (t - 5^(3))^(1) * (t - 5^(1)) - f, y * (t - 5^(1))^(2) - g )
eliminate(I, t)
152587890625*x^4*y^2-610351562500*x^3*y^2-35728465575000000000*x^3*y
      +915527343750*x^2*y^2-7444919166800000000*x^2*y-610351562500*x*y^2-
      300066291800000000*x*y+152587890625*y^2-893802369030533349376*x



mons = DynamicPolynomials.monomials((1+x)^10*(1+y)^8)
f = dot(exp.(10*rand(length(mons))),mons)

f = 421875000*x^2*y+54474750*x*y+29791*y^2-476379541*x

mons = DynamicPolynomials.monomials(f)
exps = exponents.(DynamicPolynomials.terms(f))
coeffs = DynamicPolynomials.coefficients.(f)
vtcs = hcat(exps...)'
#weights = rationalize.(log10.(abs.(coeffs)))
QQ5 = FlintPadicField(5,8)

weights = valuation.(QQ5.(coeffs)).*(-1)
p, maxcells, vertsoftrop, upperfacets = plot2Dsubdivision(vtcs,weights)
plot(p)

p = plotTropicalCurve(vertsoftrop,vtcs,upperfacets,maxcells)

plotUpperHull(vtcs, weights)

A = amoeba(f);
plot(A)
amoeba_spine_plot = plot!(spine(f))

