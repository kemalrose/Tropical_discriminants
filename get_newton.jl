
using Oscar

function get_Polytope(V_0, getVtx)
    println("get_Polytope")
    V_old, F_old = [], []
    V_new  = copy(V_0)
    F_new = []
    P = convex_hull(hcat(V_0...)')
    
    while size(V_old) != size(V_new)
        V_old = copy(V_new)
        F_old = copy(F_new)
        #println("in while loop")
        #println(V_old,F_old,getVtx)
        V_new, F_new, P = update(V_old, F_old, getVtx)
    end
    return copy(V_new), F_new, P
end

function update(V_old, F_old, getVtx)
    println("-------- UPDATING ...")
    V_new = copy(V_old)
    P = convex_hull(hcat(V_old...)')
    F_new  = [Vector{Int64}(f.a)  for f in facets(P)]
    keep_facets = []
    println(" setdiff: $(setdiff(F_new,F_old))")
    println(" F_new: $(F_new)")
    println(" F_old: $(F_old)")


    for a in setdiff(F_new, F_old)
        vert = getVtx(a)
        if !in(vert, V_new)
        #if setdiff(V_old,vert) == V_old
            println("NEW VERTEX")
            println("------ V_old : $V_old")
            println("------ vert : $vert")
            V_new = push!(V_new, vert)
        else
            println("NO NEW VERTEX -> FOUND FACET")
            println("------ V_old : $V_old")
            println("------ vert : $vert")
            index = findfirst(ℓ -> ℓ == a, F_new)
            keep_facets = push!(keep_facets, index)
            println("keep_facets: $keep_facets")
        end
    end
    #println(V_new)
    println("List of facets : $([copy(F_old); F_new[keep_facets]])")
    println("------- DONE UPDATING ")

    copy(V_new), [copy(F_old); F_new[keep_facets]], P
end
