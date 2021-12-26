
using LinearAlgebra, HomotopyContinuation, Oscar

A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]

A = [1 1 1 1 1 1;
0 1 2 0 1 0;
0 0 0 1 1 2]

A =
[
3 2 2 1 1 1 0 0 0 0;
0 1 0 2 1 0 3 2 1 0;
1 1 1 1 1 1 1 1 1 1]

A = [1 1 1 1 1 1  1 1 1 1 1 1;
      0 0 0 0 0 0  1 1 1 1 1 1;
      0 0 0 1 1 2 0 0 0 1 1 2;
      0 1 2 0 1 0 0 1 2 0 1 0]


data = data_from_matrix(A)

nsamples = 20
n = data.n
for i = 1:nsamples
    w = rand(-10000:10000,n)

    is_contained, is_in_face = cone_containments(w, data)
    non_genericity = sum(is_in_face.==1)
    println("is not generic: $non_genericity")
    #monomial = getVertex(w, data)
    #monomials = unique!(push!(monomials,monomial))
end


verify(disc, mons, coeff, A, B; strategy = "interpolation")






function interpolate_discr(A)
    data = data_from_matrix(A)
    println("Auxillary data constructed")
    v_0, vtcs, fcts, Pol = newton_pol(data)
    println("Newton polytope constructed")
    L_pts_proj = lattice_points(Pol)
    mons = [ U_inv * [v ; zeros(d)] + v_0 for v in Vector{Int64}.(L_pts_proj)]
    mons = [ [m.num for m in mon] for mon in mons]
    println("Lattice points computed")
    pts = [Horn_param(data,randn(ComplexF64, n-d),randn(ComplexF64,d)) for j = 1:length(mons)*2]
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

disc, B, intcoeff, coeff, mons, err_tol = interpolate_discr(A)


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
    (intcoeff'monomials)[1], B, intcoeff, coeff, mons, err_tol
end

#println(V_new)
#rand_verts  = sampleRandom(Γ, dets, RR, Fσ, 1000)
#P = convex_hull(  hcat(rand_verts...)'  )
#l_points  = lattice_points(P)
#verts = vertices(P)
#
#rand_verts_proj = [Π*(v-rand_verts[1]) for v in rand_verts]
##rand_verts_proj = [Π*(v-v_0) for v in rand_verts]
#P_proj = convex_hull( (hcat(rand_verts_proj...)'))
#l_points_proj  = lattice_points(P_proj)
#verts_proj = vertices(P_proj)
