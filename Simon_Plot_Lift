
function orderVtcs(vtcs)
    meanvert = mean(vtcs,dims=1)
    diffvecs = vtcs.-meanvert
    println(diffvecs)
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
    for i = 1:size(R,1)
        if R[i][end] < 0
            w = [dot(R[i],lifted_points[j,:]) for j = 1:size(lifted_points,1)]
            vertsincell = findall(j->w[j]==minimum(w),1:length(w))
            push!(maxcells,vertsincell)
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
    return p, maxcells
end



function plotUpperHull(vtcs,weights; hullcolor = "black", groundcolor = "blue")
    p, maxcells = plot2Dsubdivision(vtcs,weights)
    i = 1
    v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
    v = convert.(Int64,hcat(collect(vertices(convex_hull(v)))...)')
    v = orderVtcs(v)
    q = plot([v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = mycolor,legend = false)
    plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = groundcolor)
    for i = 1:length(maxcells)
        v = [vtcs[maxcells[i],:] weights[maxcells[i]]]
        v = convert.(Int64,hcat(collect(vertices(convex_hull(v)))...)')
        v = orderVtcs(v)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],[v[:,3];v[1,3]],color = mycolor,legend = false)
        plot!(q,[v[:,1]; v[1,1]],[v[:,2];v[1,2]],zeros(size([v[:,1]; v[1,1]])),color = groundcolor)
    end
    return q
end

#=


using Plots, Statistics, Oscar, DynamicPolynomials
@polyvar x y
f = exp(1)*x^2*y + exp(3)*y^2 + exp(8)*x^2*y^3 + exp(1)*y^4 + exp(1)*x^4*y^4
exps = exponents.(DynamicPolynomials.terms(f))
coeffs = DynamicPolynomials.coefficients.(f)
vtcs = hcat(exps...)'
weights = rationalize.(log.(coeffs))


p, maxcells = plot2Dsubdivision(vtcs,weights)
plot(p)
plotUpperHull(vtcs,weights)
=#