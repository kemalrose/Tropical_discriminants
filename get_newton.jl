
using Oscar

function get_Polytope(V_0, getVtx)
    V_old, F_old   = [], []
    V_new  = V_0
    F_new = []

    while size(V_old) != size(V_new)
        V_old = V_new
        F_old = F_new

        V_new, F_new = update(V_old, F_old, getVtx)
    end

    return V_new
end

function update(V_old, F_old, getVtx)
    V_new = V_old
    P = convex_hull(hcat(V_old...)')
    F_new  = [Vector{Int64}(f.a)  for f in facets(P)]
    keep_facets = []

    for a in setdiff(F_new, F_old)
        vert = getVtx(a)
        if ! in(vert, V_old)
            push!(V_new, vert)
        else
            index = findfirst(ℓ -> ℓ == a, F_new)
            push!(keep_facets, index)
        end
    end
    V_new, [F_old; F_new[keep_facets]]
end


V_0 = [[1 ;1], [1; 0], [0; 1]]
function getV(u)
    P = convex_hull([hcat(V_0...)'; [-1 0; 0 -1]])
    verts  = [Vector{Int64}(v) for v in vertices(P)]
    inner_prod  = [ (u' * v) for v in verts ]
    verts[argmax(inner_prod)]
end


V_old, F_old   = [], []
V_new  = V_0
V_old = V_new
V_new, F_new = update(V_old, F_old, getV)



get_Polytope(V_0, getV)
