
using Oscar#master

using Plots

function trop_curve(θ, m, n, p)
    k = length(θ)
    A = zero(MatrixSpace(QQ, k+1, 2))
    d = maximum([sum(m), sum(n)])

    for i in 1:k
        A[i, :] = [1 -θ[i]]
    end
    A[end, :] = [0 1]


    B = zero(MatrixSpace(QQ, 2, k+1))
    B = zeros(fmpq, 2, k+1)
    for i in 1:k
        B[1, i] = m[i]
        B[2, i] = n[i]
    end
    B[1, end] = d - sum(m) 
    B[2, end] = d - sum(n) 

    cycle  = trop_lin_space(A, p)
    verts = cycle.PROJECTIVE_VERTICES
    verts = zero(MatrixSpace(QQ, size(verts, 1), size(verts, 2)))
    verts = zeros(fmpq, size(verts,1), size(verts,2))
    for i in 1:size(verts, 1) 
        for j in 1:size(verts, 2)
          verts[i, j] = Rational{Int64}(cycle.PROJECTIVE_VERTICES[i, j])
        end
    end

    verts_projected = hcat(verts[:,1], verts[:, 2:end] * transpose(B))
    max_pols = [ findall(cycle.MAXIMAL_POLYTOPES[i, :]) for i in 1:size(cycle.MAXIMAL_POLYTOPES, 1)]
    verts_projected, max_pols
end



function trop_lin_space(A::fmpq_mat, p)
    QQp = FlintPadicField(p,8)
    matroid = Polymake.matroid.Matroid(VECTORS = A)
    A_val = matrix(QQp,fmpq.(A))
    vals = [ valuation(det(A_val[base.+1, :])) for base in matroid.BASES]
#    vals .= 0
    Polymake.Shell.vals = Polymake.Vector{Polymake.Rational}(vals)
    matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals, VECTORS = A)
    cycle = Polymake.tropical.linear_space(matroid_val)
    
    cycle
end



function plot_curve(verts, max_pols)
    plt = plot([])
    for i in 1:length(max_pols)
        pol = max_pols[i]
        v1 = verts[pol[1],:]
        v2 = verts[pol[2],:]
        if v1[1] + v2[1] == 2
            pt1, pt2 = v1[2:end], v2[2:end] 
        elseif v1[1] == 1
            pt1, pt2 = v1[2:end], v1[2:end] + 5 * v2[2:end] 
        elseif v2[1] == 1
            pt1, pt2 = v2[2:end], v2[2:end] + 5 * v1[2:end] 
        end
        pt1 = convert(Vector{Rational{Int64}}, pt1)
        pt2 = convert(Vector{Rational{Int64}}, pt2)
    
        plot!(plt, [pt1[1], pt2[1]], [pt1[2], pt2[2]])
    end
    plt
end



