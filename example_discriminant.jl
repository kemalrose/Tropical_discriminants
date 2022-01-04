
using Arpack, Oscar, HomotopyContinuation, LinearAlgebra
import Pkg; Pkg.add("MultiFloats")
using MultiFloats
using GenericSVD
MultiFloats.use_bigfloat_transcendentals()


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

A =
[1  1  1  1  0  0  0  0;
 0  0  0  0  1  1  1  1;
 2  3  1  0  1  0  2  0;
 2  0  0  1  0  2  2  1;
]

A =
[
  1   1   1   1   0   0   0   0 ;
  0   0   0   0   1   1   1   1;
  2   3   5   7  11 13 17 19;
19  17 13 11   7   5   3   2;
]

A = A[:, 2:end]

data = data_from_matrix(A)
n, d = data.n, data.d
v_0, vtcs, fcts, Pol = newton_pol(data)
Pol_ambient = convex_hull(hcat([data.Π_rinv * v + v_0 for v in vtcs] ...)')
Pol = Pol_ambient

Lpts = lattice_points(Pol)
fvec = f_vector(Pol)
Flist = [faces(Pol, i) for i in 0:Oscar.dim(Pol)]
LptList  = [[lattice_points(F) for F in list] for list in Flist]
IntLptList = [[interior_lattice_points(F) for F in list] for list in Flist]





A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]
#A =[3 2 2 1 1 1 0 0 0 0;0 1 0 2 1 0 3 2 1 0;1 1 1 1 1 1 1 1 1 1]
Δ, intcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x2);
verify(Δ, mons, coeff, A, B; strategy = "substitution")




disc, intcoeff, coeff, mons, err_tol = interpolate_discr(A; interpolation_method = "SVD")

verify(disc, mons, coeff, data; strategy = "interpolation")


[Horn_param(A,B,rand(-10:10, n-d)//10,rand(-10:10,d)//10) for j = 1:length(mons)*1.2]

n, d = data.n, data.d
@var λ[1:d], v[1:n-d]
φλ = [prod(λ.^A[:,i]) for i in 1:n]
Φ = φλ.*(B*v)
expand(subs(disc, variables(disc) => Φ)) == 0


V_ambient = sampleRandom(data, 2(n - d))
while Oscar.dim(convex_hull(hcat(V_ambient...)')) != n - d
    V_ambient = unique!([V_ambient; sampleRandom(data,(n - d))])
end
Pol = convex_hull(hcat(V_ambient...)')
#nsamples = 20


V = sampleRandom(data, 7000)
Pol = convex_hull(hcat(V...)')
Oscar.dim(Pol)



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
