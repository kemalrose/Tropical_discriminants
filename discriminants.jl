

using Oscar, LinearAlgebra, HomotopyContinuation


function rk_one_matrices(n::Int)
    A = zero_matrix(ZZ,2n,n^2)
    column_index  = 0
    for i in 1:n
        for j in 1:n
            column_index += 1
            A[i, column_index] += 1
            A[n + j, column_index] += 1
        end
    end
    A
end

function ishomogeneous(A)
    (n, m) = size(A)
    T = eltype(A)
    v = ones(Int64, m)
    rank(A) == rank([A; transpose(v)])
end

function tropicalized_discriminant(A)
    #if ! ishomogeneous(A)
    #    error("Not a homogeneous toric variety.")
    #end

    n, m = size(A)
    S = MatrixSpace(ZZ, n, m)
    B = nullspace(S(A))[2]
    A = Matrix{Int64}(A)
    B = Matrix{Int64}(B)

    matroid = Polymake.matroid.Matroid(VECTORS = B)
    ΣB = Polymake.tropical.matroid_fan{min}(matroid)

    max_cones = []
    raymatrix = Matrix{Int64}(ΣB.PROJECTIVE_VERTICES)[1:end-1,2:end]
    raymatrix = vcat(raymatrix, A, -A)
    n_rays = size(raymatrix, 1)

    incidence_max_cones = Matrix{Int64}(ΣB.MAXIMAL_POLYTOPES)[:, 1:end-1]
    n_cones = size(incidence_max_cones,1)
    incidence_max_cones = hcat(incidence_max_cones, ones(Int64, n_cones, 2*n))

    for ind in 1:size(incidence_max_cones, 1)
        ray_list = filter(λ -> Bool(incidence_max_cones[ind, λ]), 1:n_rays)
        cone = positive_hull(raymatrix[ray_list, :])
        push!(max_cones, Polyhedron(facets(cone),linear_span(cone)))
        #push!(max_cones,cone)
    end
    max_cones
end





function ray(ω, i)
    n = length(ω)
    A = Matrix{Int}(vcat( -identity_matrix(ZZ, n)[ filter(λ -> λ != i, 1:end ), :], identity_matrix(ZZ, n)))
    b = vcat(-ω[filter(λ -> λ != i, 1:n)], ω)
    Polyhedron(A, b)
end

tropical_multiplicity(A, cone)
    1
end

S = MatrixSpace(ZZ, 2, 3)
#A = [1 0 2 1; 0 1 -1 0]
A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]


ω = [12, 43, 22, 1, 6, 2]
@var a[1:6]
monomial = 1
maxcones = tropicalized_discriminant(A)
for i in 1:6
    Ray = ray(ω, i)

    [intersect(Ray, cone) for cone in maxcones].>-1

end
