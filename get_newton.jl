
using Oscar

function Horn_param(A, B, v, λ)
    φλ = [prod(λ.^A[:,i]) for i in 1:n]
    φλ.*(B*v)
end

function get_Vdm(pts, mons; normalize = true)
    T = eltype(pts[1])
    V = zeros(T, (length(pts), length(mons)))
    for i in 1:length(pts)
        V[i, :] = [prod(pts[i].^mon) for mon in mons]
        if normalize
            V[i, :] = V[i, :]/norm(V[i, :])
        end
    end
    V
end

for mon in mons
    println(mon)
    x = prod(pts[25].^mon)
end

function interpolate_discr(A)
    d, n = size(A)
    S = MatrixSpace(ZZ, d, n)
    B = nullspace(S(A))[2]
    S,U,V = snf_with_transform(B)
    Π = Matrix{Int64}(U[1:n-d, :])
    A = Matrix{Int64}(A)
    B = Matrix{Int64}(B)
    U_inv  = Matrix{Int64}(inv(U))


    v_0, vtcs, fcts, Pol = newton_pol(A,B,Π)
    L_pts_proj = lattice_points(Pol)
    mons = [ U_inv * [v ; zeros(d)] + v_0 for v in Vector{Int64}.(L_pts_proj)]
    mons = [ [m.num for m in mon] for mon in mons]

    pts = [Horn_param(A,B,randn(ComplexF64, n-d),randn(ComplexF64,d)) for j = 1:length(mons)*2]

    V = get_Vdm(pts, mons)

    U, sings, VV = svd(V)
    coeff = VV[:, end]

    max_val  = maximum(abs.(coeff))
    coeff = coeff./ coeff[findfirst(k-> abs(coeff[k]) == max_val, 1:length(coeff))]
    coeff = real.(coeff)


    err_tol = sings[end]/sings[end-1] * 100
    coeff = rationalize.(coeff, tol = err_tol)

    @var a[1:n]
    monomials = [prod(a.^mon) for mon in mons]
    (coeff'monomials)[1], B, coeff, mons, err_tol
end



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
