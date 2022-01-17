
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


#BBernd M
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
Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x2, sample_method = "Horn", redundancy_factor = 1.5);
bool = verify(Δ, mons, coeff, data; strategy = "substitution")

coeff_fin, mons_fin, F = interpolate_discr_finitely(A, redundancy_factor = 1.5)


function compare_results(ratcoeff, coeff_fin, F)
    denom = lcm([c.den for c in ratcoeff])
    intcoeff = Int64.(ratcoeff.*denom)
    rats = [ F(coeff) ÷ F(intcoeff[1]) for coeff in intcoeff]
    rats2 = [ coeff ÷ coeff_fin[1] for coeff in coeff_fin]
    sum((rats-rats2).!=0) == 0
end


function reparameterise(A)
    rk = rank(A)
    d, n = size(A)
    S,T,U = AbstractAlgebra.snf_with_transform(matrix(ZZ, A))
    TT = Matrix{Int64}(Matrix(inv(T)))
    hcat(TT[:, 1:rk], zeros(Int64, (d, n-rk)) ) * Matrix{Int64}(Matrix(inv(U)))
end





 A =
     [ 1 1 1 1 0 0 0 0;
     0 0 0 0 1 1 1 1;
     2 3 5 7 11 13 17 19;
     19 17 13 11 7 5 3 2]

A = A[:, 1:end.!=4]

for i in 1:8
    println("________________________")
    println("i = $i")
     A_col = A[:, 1:end.!=i]
     A_col = reparameterise(A_col)
     println("Index = ", index(A_col))
     v_0, vtcs, fcts, Pol = newton_pol(A_col)
     n_mons  = length(lattice_points(Pol))
     println("Nr of monomials = $n_mons")
end
comput_result  = []
for i in 1:8
    println("________________________")
    println("i = $i")
    A_col = reparameterise(A[:, 1:end.!=i])
    println("Compute the Newton Polytope!")
    v_0, vtcs, fcts, Pol = newton_pol(A_col)
    n_mons  = length(lattice_points(Pol))
    println("There are $n_mons monomials")
    println("Now compute the discriminant!")
    Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A_col; interpolation_method = "SVD", T = Float64x8, sample_method = "Horn", redundancy_factor = 1.2)
end

n_mons = [ size(result[3], 1) for result in comput_result  ]


println("Delete column i = $i")
println("Index of AA = ", index(AA))
println("Theren are $n_mons monomials.")


filename = "/Users/kemal/Desktop/mons.txt"
filename = "/Users/kemal/Desktop/Horn.txt"

# writing to files is very similar:
f = open(filename, "w")
# both print and println can be used as usual but with f as their first arugment

println(f, "A = $A")
println(f, "B = $B")
println(f, "λ^A = [prod(λ.^data.A[:,i]) for i in 1:7]")
println(f, "(v, λ) ⟼ λ^A .* B*v")

for mon in mons
    println(f, mon)
end
close(f)


A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]

A = A[:, 1:end.!=4]


Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x8, sample_method = "Horn", redundancy_factor = 5)
bool = verify(Δ, mons, coeff, data; strategy = "substitution")

coeff, mons, F = interpolate_discr_finitely(A, p = 373)

vec_int = [12, 64, 2, 46, 3, 17, 62]
vec_F = F.(vec_int)

vec_new = F(rand(1:373-1)).*vec_F




using Plots


function xx(a)
     maximum([5+2*a;5+a;2])
end
function yy(a)
     maximum([-2+2*a;2+a;1])
end


as = -10:0.01:10


plot(xx.(as),yy.(as))


function ϕ(a)
     [10. ^(5+2*a).+10. ^(5+a).+10^2; (10).^(-2. +2*a).+10. ^(2. +a).+10^1]
end
function ψ(a)
     r = 10*randn(ComplexF64)
     [10^5*r^(2*a)+10^5*r^a+10^2;10^(-2)*r^(2*a)+10^2*r^a+10^1]
end

vals2 = log10.(abs.((transpose(hcat(ψ.(as)...)))));


plot!(vals2[:,1],vals2[:,2],seriestype = :scatter)


[log10(1) + vals2[14,:]'*[2,0]
log10(20000000) + vals2[14,:]'*[1,1]
log10(100000000000000) + vals2[14,:]'*[0,2]
log10(9998980000200) + vals2[14,:]'*[1,0]
log10(9799002000000000) + vals2[14,:]'*[0,1]
log10(8899101999990000)]

x^2  - 20000000*x*y + 100000000000000*y^2  - c*x +
 9799002000000000*y - 8899101999990000

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
