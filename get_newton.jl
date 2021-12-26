
using Oscar

function Horn_param(data, v, λ)
    φλ = [prod(λ.^data.A[:,i]) for i in 1:n]
    φλ.*(data.B*v)
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


function newton_pol(data::Aux_data)
    n, d = data.n, data.d
    V_ambient = sampleRandom(data, 2(n - d))
    while Oscar.dim(convex_hull(hcat(V_ambient...)')) != n - d
        V_ambient = unique!([V_ambient; sampleRandom(data,(n - d))])
    end
    #V_ambient = [ [vi.num for vi in v] for v in V_ambient ]
    v_0 = V_ambient[1]
    V_0 = [data.Π*(v - v_0) for v in V_ambient]

    function getV(w)
        w_new = data.Π'*w
        w_new  = w_new * 10000 + rand(-100:100, size(w_new))
        vtx = getVertex(-w_new, data)
        monomial = data.Π * (vtx - v_0)
        monomial
    end
    vtcs, fcts, Pol = get_Polytope(V_0, getV);

    v_0, vtcs, fcts, Pol
end



function get_Polytope(V_0, getVtx)
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
    copy(V_new), [copy(F_old); F_new[keep_facets]], P
end


function verify(disc, mons, coeff, data; strategy = "substitution")
    A, B, d, n = data.A, data.B, data.d, data.n

    if strategy == "substitution"
        @var λ[1:d], v[1:n-d]
        φλ = [prod(λ.^A[:,i]) for i in 1:n]
        Φ = φλ.*(B*v)
        return expand(subs(disc, variables(disc) => Φ)) == 0
    elseif strategy == "interpolation"
        pts = [Horn_param(A,B,rand(-10:10, n-d)//10,rand(-10:10,d)//10) for j = 1:length(mons)*1.2]
        pts = Vector{Rational{BigInt}}.(pts)
        V = get_Vdm(pts, mons, normalize = false)
        S = MatrixSpace(QQ, size(V,1), size(V,2))
        if rank(S(V)) == size(V, 2) - 1
            return norm(V * coeff) == 0
        else
            print("Vdm is not of full rank!")
        end
    else
        print("Choose input strategy to either be substitution or interpolation")
    end
end





function interpolate_discr(A; interpolation_method = "SVD")
    data = data_from_matrix(A)
    n, d = data.n, data.d

    println("1. Compute the Newton polytope P of the discriminant")
    v_0, vtcs, fcts, Pol = newton_pol(data)

    println("2. Find its lattice points")
    L_pts_proj = lattice_points(Pol)

    println("   P has $(length(L_pts_proj)) lattice points")
    mons = [ data.Π_rinv * v + v_0 for v in Vector{Int64}.(L_pts_proj)]
    mons = [ [m.num for m in mon] for mon in mons]

    println("3. Construct the interpolation problem")
    pts = [Horn_param(data,randn(ComplexF64, n-d),randn(ComplexF64,d)) for j = 1:convert(Int64,round(length(mons)*1.2))]
    V = get_Vdm(pts, mons)
    println("   Constructed a Vandermonde matrix of size $(size(V))")

    println("4. Find floating point coefficients...")
    if interpolation_method == "SVD"
        println("   ... using SVD")
        U, sings, VV = svd(V)
        coeff = VV[:, end]
        err_tol = sings[end]/sings[end-1] * 100
    elseif interpolation_method == "eigs"
        println("   ... using eigs")
        VtV = V' * V
        @time EVobj = eigs(VtV; nev = 2, which = :SM)
        coeff = EVobj[2][:, 1]
        println(EVobj[1])
        err_tol = abs(EVobj[1][1]/EVobj[1][2]) * 100
    else
        println("Choose interpolation_method to be either SVD or eigs.")
    end

    println("5. Rationalizing with tolerance $(err_tol)")
    nonzinds = findall(ℓ->abs(ℓ)>err_tol, coeff)
    coeff = coeff/coeff[findfirst(ℓ->abs(ℓ) == minimum(abs.(coeff[nonzinds])),coeff)]
    intcoeff = convert.(Int64,round.(real.(coeff)))
    #coeff = rationalize.(coeff, tol = err_tol)

    @var a[1:n]
    monomials = [prod(a.^mon) for mon in mons]
    (intcoeff'monomials)[1], intcoeff, coeff, mons, err_tol
end
