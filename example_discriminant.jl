using Oscar
using LinearAlgebra

A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]
n, m = size(A)
S = MatrixSpace(ZZ, n, m)
B = nullspace(S(A))[2]
S,U,V = snf_with_transform(B)
Π = Matrix(U[1:m-n, :])
U_inv  = Matrix(inv(U))

A = Matrix{Int64}(A)
B = Matrix{Int64}(B)





matroid = Polymake.matroid.Matroid(VECTORS = B)
ΣB = Polymake.tropical.matroid_fan{min}(matroid)
#ΣB = Polymake.tropical.matroid_fan_from_flats{min}(matroid)
Γ, dets, RR, Fσ = getConeData(A, ΣB)




function getV(w)
    #println(w)
    w_new = U_inv * [w; zeros(Int,n)]
    vtx, flag = getVertex2(w_new, Γ, dets, RR, Fσ)
    println("FLAG : $flag")
    Π * (vtx - v_0)
end
w_new  = [1 1 1 1 1 1]'

V_0 = sampleRandom(Γ, dets, RR, Fσ, 2(m - n))
while dim(convex_hull(hcat(V_0...)')) != m - n
    V_0 = sampleRandom(Γ, dets, RR, Fσ, 2(m - n))
end
v_0 = V_0[1]
V_0 = [Π*(v - v_0) for v in V_0]
vtcs = get_Polytope(V_0, getV);
length(vtcs)
length(unique!(vtcs))
P = convex_hull(hcat(vtcs...)')


rand_verts  = sampleRandom(Γ, dets, RR, Fσ, 1000)
P = convex_hull(  hcat(rand_verts...)'  )
l_points  = lattice_points(P)
verts = vertices(P)

rand_verts_proj = [Π*(v-rand_verts[1]) for v in rand_verts]
rand_verts_proj = [Π*(v-v_0) for v in rand_verts]
P_proj = convex_hull( (hcat(rand_verts_proj...)'))
l_points_proj  = lattice_points(P_proj)
verts_proj = vertices(P_proj)










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
