
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

A =      [ 1 1 1 1 0 0 0 0;
       0 0 0 0 1 1 1 1;
       2 3 5 7 11 13 17 19;
      19 17 13 11 7 5 3 2]



data = data_from_matrix(A)

disc, intcoeff, coeff, mons, err_tol = interpolate_discr(A)

verify(disc, mons, coeff, data; strategy = "interpolation")

#nsamples = 20
#n = data.n
#for i = 1:nsamples
#    w = rand(-10000:10000,n)
#
#    is_contained, is_in_face = cone_containments(w, data)
#    non_genericity = sum(is_in_face.==1)
#    println("is not generic: $non_genericity")
#    #monomial = getVertex(w, data)
#    #monomials = unique!(push!(monomials,monomial))
#end





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
