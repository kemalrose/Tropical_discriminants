
using Oscar

function newton_pol(A, B, Π)
    d, n = size(A)

    matroid = Polymake.matroid.Matroid(VECTORS = B)
    ΣB = Polymake.tropical.matroid_fan{min}(matroid)
    #ΣB = Polymake.tropical.matroid_fan_from_flats{min}(matroid)
    Γ, dets, RR, Fσ = getConeData(A, ΣB)


    V_0 = sampleRandom(Γ, dets, RR, Fσ, 2(n - d))
    while dim(convex_hull(hcat(V_0...)')) != n - d
        V_0 = sampleRandom(Γ, dets, RR, Fσ, 4(n - d))
    end
    #V_0 = [ [vi.num for vi in v] for v in V_0 ]
    v_0 = V_0[1]
    V_0 = [Π*(v - v_0) for v in V_0]

    function getV(w)
        w_new = Π'*w
            w_new  = w_new * 10000 + rand(-100:100, size(w_new))
        vtx, flag = getVertex2(-w_new, Γ, dets, RR, Fσ)
        Π * (vtx - v_0)
    end
    vtcs, fcts, Pol = get_Polytope(V_0, getV);

    v_0, vtcs, fcts, Pol
end



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
    V_new = copy(V_old)
    P = convex_hull(hcat(V_old...)')
    F_new  = [Vector{Int64}(f.a)  for f in facets(P)]
    keep_facets = []



    for a in setdiff(F_new, F_old)
        vert = getVtx(a)
        if !in(vert, V_new)
        #if setdiff(V_old,vert) == V_old
            V_new = push!(V_new, vert)
        else

            index = findfirst(ℓ -> ℓ == a, F_new)
            keep_facets = push!(keep_facets, index)
        end
    end
    #println(V_new)

    copy(V_new), [copy(F_old); F_new[keep_facets]], P
end
