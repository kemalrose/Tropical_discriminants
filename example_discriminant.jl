using Oscar
using LinearAlgebra

A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]

A = [1 1 1 1 1 1;
0 1 2 0 1 0;
0 0 0 1 1 2]

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
mons = [ U_inv * [v ; zeros(n-d)] + v_0 for v in Vector{Int64}.(L_pts_proj)]
mons = [ [m.num for m in mon] for mon in mons]

pts = [Horn_param(A,B,randn(F64,n-d),randn(F64,d)) for j = 1:50]
pts_proj = [[prod(p.^Π[i,:]) for i in 1:n-d] for p in pts]
L_points = Vector{Int64}.(L_points)

V = get_Vdm(pts, mons)
coeff = nullspace(V)
svd(V)



function Horn_param(A, B, v, λ)
    φλ = [prod(λ.^A[:,i]) for i in 1:n]
    φλ.*(B*v)
end

function get_Vdm(pts, mons)
    V = zeros(Complex{Float64}, (length(pts), length(mons)))
    for i in 1:length(pts)
        V[i, :] = [prod(pts[i].^mon) for mon in mons]
        V[i, :] = V[i, :]/norm(V[i, :])
    end
    V
end

function interpolate_discr(A)
    B =
    P =
    mons = lattice_points(P)

    n_pts = length(mons) * 2

    V
end




x = pts[1]
det([2*x[1] x[2] x[4]; x[2] 2*x[3] x[5]; x[4] x[5] 2*x[6]])

for x in pts
    println(abs(det([2*x[1] x[2] x[4]; x[2] 2*x[3] x[5]; x[4] x[5] 2*x[6]])
))
end





function getV(w)
    #println(w)
    w_new = Π'*w
    w_new  = w_new * 1000 + rand(1:10, size(w_new))
    vtx, flag = getVertex2(w_new, Γ, dets, RR, Fσ)
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
vertices(P)





rand_verts  = sampleRandom(Γ, dets, RR, Fσ, 1000)
P = convex_hull(  hcat(rand_verts...)'  )
l_points  = lattice_points(P)
verts = vertices(P)

rand_verts_proj = [Π*(v-rand_verts[1]) for v in rand_verts]
#rand_verts_proj = [Π*(v-v_0) for v in rand_verts]
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



function getV_proj(u)
    verts  = [Vector{Int64}(v) for v in vertices(P_proj)]
    inner_prod  = [ (u' * v) for v in verts ]
    verts[argmax(inner_prod)]
end

function getV(u)
    verts  = [Vector{Int64}(v) for v in vertices(P)]
    inner_prod  = [ (u' * v) for v in verts ]
    verts[argmax(inner_prod)]
end

function adapted_getV(u)
     u_new = Π' * u
     #u_new  = 100*u_new + rand(1:10, size(u_new))
     vert = getV(u_new)
     Π * (vert - rand_verts[1])
end



function getV(w)
    w_new = Π'*w
    w_new  = w_new * 10000 + rand(-100:100, size(w_new))
    vtx, flag = getVertex2(-w_new, Γ, dets, RR, Fσ)
    println("w : $w")
    println("w_new : $w_new")
    println("flag : $flag")
    Π * (vtx - v_0)
end

v_0 = rand_verts[1]
V_0 = rand_verts_proj[1:5]
V_old, F_old   = [], []
V_new  = V_0
V_old = V_new
V_new, F_new = update(V_old, F_old, getV_proj)

verts, fcts, P = get_Polytope(V_0, getV)

L_points  = lattice_points(P)
