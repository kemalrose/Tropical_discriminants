



function tropicalize(f)
    mons = DynamicPolynomials.monomials(f)
    exps = exponents.(DynamicPolynomials.terms(f))
    coeffs = DynamicPolynomials.coefficients.(f)
    τf = vcat([ [log10(abs(coeffs[i])) exps[i][1]] for i in 1:length(mons)]...)
    τf
end


function collect_roots(τf)
    trop_roots = []
    mults = []
    for i in 1:size(τf, 1)
        for j in i+1:size(τf, 1)
            l = τf[i, :] - τf[j, :]
            root = -l[1]/l[2]
            max = maximum( τf * [1;root] )
            #Have to check that only two functions "win"
            if abs(max - (τf*[1;root])[i]) < 1e-5
                push!(trop_roots, root)
                push!(mults, Int64(abs(l[2])))
            end            
        end
    end
    #Need that all gathered roots are unique!
    σ = sortperm(trop_roots)
    trop_roots[σ], mults[σ]
end





function appr_trop_curve(f, g)  
    τf, τg = tropicalize(f), tropicalize(g)
    roots_τf, mults_f = collect_roots(τf)
    roots_τg, mults_g = collect_roots(τg)
    all_roots = append!([], roots_τf, roots_τg )
    all_mults = append!([], mults_f, mults_g)
    #Assume all roots of τf and τg to be distinct
    σ = sortperm(all_roots)
    all_roots = all_roots[σ]
    all_mults = all_mults[σ]

    pts_between = append!([all_roots[1] - 1], [ (all_roots[i] + all_roots[i+1])/2. for i in 1:length(all_roots)-1], [all_roots[end]+1])
    
    xmin = minimum([ maximum(τf * [1; a]) for a in all_roots]) - 2 
    ymin = minimum([ maximum(τg * [1; a]) for a in all_roots]) - 2
    
    τf_vals = [ maximum(τf * [1; a]) for a in all_roots]
    τg_vals = [ maximum(τg * [1; a]) for a in all_roots]
    plot(τf_vals, τg_vals, colour = "green", legend = false)
    
    vtcs = [ [τf_vals[i];τg_vals[i]] for i in 1:length(τf_vals) ]
        

    slopes = []
    for pt in pts_between
        index_f  = argmax(τf * [1; pt])
        dx = τf[index_f,2]

        index_g  = argmax(τg * [1; pt])
        dy = τg[index_g,2]

        push!(slopes,[dx; dy])
    end

    nverts = length(vtcs)
    edges = []
    edge_mults = []
    for i in 1:nverts-1
        push!(edges, [i, i+1])
        push!(edge_mults, 0)
    end
    
    for i in 1:nverts
        slope = slopes[i] - slopes[i+1]
        vert = vtcs[i]
        if slope[2] < 0
            push!(vtcs, [vert[1]; ymin] )
            push!(edges, [i, length(vtcs)])
            push!(edge_mults, all_mults[i])
        elseif slope[1] < 0
            push!(vtcs, [xmin; vert[2]] )
            push!(edges, [i, length(vtcs)])
            push!(edge_mults, all_mults[i])
        end
    end
    
    a = all_roots[end] + 1
    push!(vtcs, [maximum(τf * [1; a]); maximum(τg * [1; a])])
    push!(edges, [nverts, length(vtcs)])
    push!(edge_mults, 0)


    for i in 1:nverts
        v = vtcs[i]
        adj_indices = findall(e -> i in e, edges)
        mults = edge_mults[adj_indices]

        b = [0; 0]
        A = [0; 0]
        for j in 1:length(adj_indices)
            e = edges[adj_indices][j]
            if e[1] == i
                s = vtcs[e[2]] - vtcs[e[1]]
            else
                s = vtcs[e[1]] - vtcs[e[2]]
            end
            s = s/maximum(abs.(s))
            rs = rationalize.(s, tol = 1e-8)
            dens = denominator.(rs)
            rs = rs.*lcm(dens)

            if mults[j] != 0
                b += mults[j] * rs
            else
                A = rs
            end
        end
        edge_mults[adj_indices[findfirst(m->m==0, mults)]] = -Int64(A\b)
    end
    

    found = false
    index  = 1
    while found == false && index <= length(vtcs)
        edge_inds = findall(e -> index in e, edges)
        if length(edge_inds) == 2
            vts = setdiff(vcat(edges[edge_inds]...), index)
            found = true
            deleteat!(edges, edge_inds)
            push!(edge_mults, edge_mults[edge_inds[1]])
            deleteat!(edge_mults, edge_inds)
            push!(edges, [vts[1], vts[2]])
            deleteat!(vtcs, index)
            for e in edges
                e[e.>index] = e[e.>index].-1
            end
        end
        index += 1
    end


    xvals = [v[1] for v in vtcs]
    yvals = [v[2] for v in vtcs]
    p =   plot(xvals[edges[1]], yvals[edges[1]], legend = false, colour = "red")
    for e in edges[2:end]
        plot!(p, xvals[e], yvals[e], colour = "red")
    end


    p, vtcs, edges, edge_mults
end




function upper_hull(f, g; hullcolor = "red", groundcolor = "blue", rangemin = 5, rangemax = 10)
    τf, τg = tropicalize(f), tropicalize(g)
    roots_τf, mults_f = collect_roots(τf)
    roots_τg, mults_g = collect_roots(τg)

    y_coords = [ maximum(τg * [1; a]) for a in roots_τf]
    x_coords = [ maximum(τf * [1; a]) for a in roots_τg]

    exps_vertical = [[0 0]]
    weights_vertical = [0.]
    for i in 1:length(roots_τf)
        n_vert = exps_vertical[end] + [0 mults_f[i] ]
        push!(exps_vertical, n_vert)
        push!(weights_vertical, weights_vertical[end] -mults_f[i]*y_coords[i])
    end
    exps_horizontal = [[0 0]]
    weights_horizontal = [0.]
    for i in 1:length(roots_τg)
        n_vert = exps_horizontal[end] + [mults_g[i] 0]
        push!(exps_horizontal, n_vert)
        push!(weights_horizontal, weights_horizontal[end] -mults_g[i]*x_coords[i])
    end

    exps = vcat(vcat(exps_vertical...), vcat(exps_horizontal[2:end]...))
    weights = [weights_vertical; weights_horizontal[2:end]]
    weights = rationalize.(weights, tol = 1e-8)
    weights =  weights .- minimum(weights) .+= rangemin

    #weights =  weights .- minimum(weights)
    #max = maximum(weights)
    #weights = weights .* (rangemax - rangemin)//(max) .+= rangemin
    plotUpperHull(exps,weights; hullcolor = hullcolor, groundcolor = groundcolor)
end

#=

@polyvar x y t

f = 10^5 * t^2 + 10^5*t + 100
g = 1/100 * t^2 + 100 * t + 10


f = sum([10^(rand()*10) *t^i for i in 0:8])
g = sum([10^(rand()*10) *t^i for i in 0:8])

p, vtcs, edges, edge_mults = appr_trop_curve(f, g)  
plot(p)
f = 10^5 * t^2 + 10^5*t + 100
g = 1/100 * t^2 + 100 * t + 10

q! = upper_hull(f, g)

p, vtcs, edges, edge_mults = appr_trop_curve(f, g)
plot(p)
=#