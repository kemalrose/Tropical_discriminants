using Oscar
using LinearAlgebra, HomotopyContinuation

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


V_0 = sampleRandom(Γ, dets, RR, Fσ, 1000)
P = convex_hull(hcat(V_0...)')
v_0, vtcs, fcts, Pol = newton_pol(A, B, Π)
disc, B, coeff, mons, err_tol = interpolate_discr(A)

verify(disc, mons, coeff, A, B; strategy = "interpolation")


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
