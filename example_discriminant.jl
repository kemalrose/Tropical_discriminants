
using Arpack, Oscar, HomotopyContinuation, LinearAlgebra
import Pkg; Pkg.add("MultiFloats")
using MultiFloats
using GenericSVD
MultiFloats.use_bigfloat_transcendentals()
#Pkg.add("https://github.com/JuliaLinearAlgebra/GenericSVD.jl.git")


A =
[1 1 1 1 1 1;
0 1 2 0 1 0;
0 0 0 1 1 2]

A =
[
3 2 2 1 1 1 0 0 0 0;
0 1 0 2 1 0 3 2 1 0;
1 1 1 1 1 1 1 1 1 1]

A =
    [1 1 1 1 1 1  1 1 1 1 1 1;
    0 0 0 0 0 0  1 1 1 1 1 1;
    0 0 0 1 1 2 0 0 0 1 1 2;
    0 1 2 0 1 0 0 1 2 0 1 0]

A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]

A =
    [1  1  1  1  0  0  0  0;
    0  0  0  0  1  1  1  1;
    2  3  1  0  1  0  2  0;
    2  0  0  1  0  2  2  1;]

A =
    [
    1   1   1   1   0   0   0   0 ;
    0   0   0   0   1   1   1   1;
    0   1   3   5  0 2 6 8;
    8  6 2 0   5   3   1   0;
    ]

A = A[:, 2:end]


A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]


AA =
[
 1 1 1 1 0 0 0;
 0 0 0 0 1 1 1;
 0 0 0 1 0 1 4;
 3 2 0 0 0 0 2;
]
A = AA



#A =[3 2 2 1 1 1 0 0 0 0;0 1 0 2 1 0 3 2 1 0;1 1 1 1 1 1 1 1 1 1]
Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x4, sample_method = "Horn", redundancy_factor = 2);
bool = verify(Δ, mons, coeff, data; strategy = "substitution")

coeff_fin, mons_fin, F = interpolate_discr_finitely(A, redundancy_factor = 1.5)


denom = lcm([c.den for c in ratcoeff])
intcoeff = Int64.(ratcoeff.*denom)

rats = [ F(coeff) ÷ F(intcoeff[1]) for coeff in intcoeff]
rats2 = [ coeff ÷ coeff_fin[1] for coeff in coeff_fin]
sum(rats-rats2)

function index(A::Matrix{Int}, S=AbstractAlgebra.snf(matrix(ZZ, A)))
    # Returns the index of the lattice defined by A if it is of full rank and 0 else.
    n, m = size(A)
    if n > m
        return 0
    end
    prod( [S[i,i] for i in 1:n] )
end



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
